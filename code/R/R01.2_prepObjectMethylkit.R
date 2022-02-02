## MethylKit object preparation
## A. Balard
## 25th of August 2021

## R/4.0.2 to run
library(methylKit)
library(readxl)
library(plyr) # for join (keep row order)
library(dendextend) # The package dendextend contains many functions for changing the appearance of a dendrogram and for comparing dendrograms.
library(ggplot2)

## Sources:
## https://www.bioconductor.org/packages/devel/bioc/vignettes/methylKit/inst/doc/methylKit.html
## https://nbis-workshop-epigenomics.readthedocs.io/en/latest/content/tutorials/methylationSeq/Seq_Tutorial.html#load-datasets
## /data/archive/archive-SBCS-EizaguirreLab/RRBS/StickPara_Broject_archive/08CompGenomes_mCextr/04RAnalyses_Methylome/methylbam_peichel...

## load custom functions
source("/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/GIT_StickParaOffsBroject/code/R/customRfunctions.R")

##### Load prepared dataset #####
dataPath="/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/Data/04BSBolt_methCall/BSBolt/MethylationCalling/Methylation_calling_splitted/formatCG4methylKit"

temp = list.files(path=dataPath,
                  pattern = ".CG4methylkit.txt",
                  full.names = T)

## Add metadata on treatments
metadata <- readxl::read_xlsx("/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/GIT_StickParaOffsBroject/data/Kostas_G2_info.xlsx")

metadata$trtG1G2_NUM <- as.numeric(as.factor(metadata$trtG1G2))

### Make methylkit object
myobj=methylKit::methRead(as.list(temp),
                          mincov=10,
                          sample.id=as.list(metadata$ID),
                          assembly="Gynogen_pchrom_assembly_all",
                          treatment=metadata$trtG1G2_NUM,
                          context="CpG")

##############################################
## We remove several samples (N=9)
## S12, S22, S110, S118 & S142 (less than 6M reads after trimming.
## some have also very weird methylation profiles, bad fastQC)
## and Fam12 (N=4, only present in parental group)
IDtoRm= c("S12", "S22", "S110", "S118", "S142",
          metadata$ID[metadata$Family %in% "Fam12"])

## create a new methylRawList object
print("Remove 9 samples")

myobj=reorganize(
    myobj,
    sample.ids=metadata$ID[!metadata$ID %in% IDtoRm],
    treatment=metadata$trtG1G2_NUM[!metadata$ID %in% IDtoRm])

###############################
## Filtering and normalising ##
###############################

## Filtering based on coverage:
# It might be useful to filter samples based on coverage. Particularly, if our samples are suffering from PCR bias it would be useful to discard bases with very high read coverage. Furthermore, we would also like to discard bases that have low read coverage, a high enough read coverage will increase the power of the statistical tests. The code below filters a methylRawList and discards bases that have coverage below 10X and also discards the bases that have more than 99.9th percentile of coverage in each sample.
print("Filter")
filtered.myobj=filterByCoverage(myobj, lo.count=10, lo.perc=NULL,
                                hi.count=NULL, hi.perc=99.9)

## normalise the coverage
print("Normalise")
normFil.myobj=normalizeCoverage(filtered.myobj)

#####################
## MERGING SAMPLES ##
#####################
##In order to do further analysis, we will need to get the bases covered in all samples. The following function will merge all samples to one object for base-pair locations that are covered in all samples. The unite() function will return a methylBase object which will be our main object for all comparative analysis. The methylBase object contains methylation information for regions/bases that are covered in all samples.

table(metadata$trtG1G2[!metadata$ID %in% IDtoRm])
#   Control  E_control  E_exposed    Exposed NE_control NE_exposed 
#        12         28         28         12         28         27 

print("Add CpG present in ALL individuals")
uniteCovALL=unite(normFil.myobj, mc.cores=8)
uniteCovALL=as(uniteCovALL,"methylBase")

## Sagonas et al. 2020. keep only those methyl-ated CpG sites observed in at least two individual fish after filtering and normalising
uniteCov2=unite(normFil.myobj, min.per.group=2L, mc.cores=8)# try with 8 cores
uniteCov2=as(uniteCov2,"methylBase")
   
## We save also DMS present in at least 6 fish
uniteCov6=unite(normFil.myobj, min.per.group=6L, mc.cores=8)# try with 8 cores
uniteCov6=as(uniteCov6,"methylBase")

## And save in RData file for later analyses
save(uniteCovALL, file = "/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/Data/05MethylKit/output/uniteCovALL.RData")
save(uniteCov2, file = "/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/Data/05MethylKit/output/uniteCov2.RData")
save(uniteCov6, file = "/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/Data/05MethylKit/output/uniteCov6.RData")

########################################################################################
## Remove reads from sex chromosome X ("Gy_chrXIX") and unmapped contigs ("Gy_chrUn") ##
########################################################################################
print("nbr CpG shared by all 137 samples:")
length(uniteCovALL$chr)

print("nbr CpG shared by at least 2 animals:")
length(uniteCov2$chr) 

print("nbr CpG shared by at least 6 animals:")
length(uniteCov6$chr) 

print("nbr CpG per chrom shared by all 137 samples:")
table(uniteCovALL$chr)

print("nbr CpG on sex chromosome of unmapped:")
nrow(uniteCovALL[uniteCovALL$chr %in% c("Gy_chrXIX", "Gy_chrUn"),])

## Keep CpG apart from sex chromosome XIX and unmapped (comprise Y chr)
uniteCovALL_woSexAndUnknowChr=uniteCovALL[!uniteCovALL$chr %in% c("Gy_chrXIX", "Gy_chrUn"),]
uniteCov2_woSexAndUnknowChr=uniteCov6[!uniteCov2$chr %in% c("Gy_chrXIX", "Gy_chrUn"),]
uniteCov6_woSexAndUnknowChr=uniteCov6[!uniteCov6$chr %in% c("Gy_chrXIX", "Gy_chrUn"),]

#######################################
## Saving point for further analyses ##
#######################################
save(uniteCovALL_woSexAndUnknowChr, file = "/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/Data/05MethylKit/output/uniteCovALL_woSexAndUnknownChr.RData")
save(uniteCov2_woSexAndUnknowChr, file = "/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/Data/05MethylKit/output/uniteCov2_woSexAndUnknownChr.RData")
save(uniteCov6_woSexAndUnknowChr, file = "/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/Data/05MethylKit/output/uniteCov6_woSexAndUnknowChr.RData")

##################### Previous tests with ALL numbers of fish 1 to 12:
## we kept for downstream analyses all CpG sites present in at least 1 to 12 individuals per group, or in all individuals:
# print("Unite and store in a list")
# mylist_uniteCov=list()
# for (i in 6:12L){ # done for 1 to 5, then 6 to 12
#     uniteCov=unite(normFil.myobj, min.per.group=i, mc.cores=8)# try with 8 cores
#     uniteCov=as(uniteCov,"methylBase")
#     name=paste0("uniteCov_", as.character(i))
#     mylist_uniteCov[[name]]=uniteCov
# }

## Add CpG present in ALL individuals
# uniteCov=unite(normFil.myobj, mc.cores=8)
# uniteCov=as(uniteCov,"methylBase")
# mylist_uniteCov[["uniteCov_ALL"]]=uniteCov

# CpGALL=length(mylist_uniteCov$uniteCov_ALL$coverage1) # 47238

# Idea: plot number of retained CpG site by nbr of individuals sharing these CpG sites. Preparing file for that:
# print("Make DF")
# CpGDF1_5=data.frame(lapply(mylist_uniteCov, function(x) length(x$coverage1)))
# CpGDF1_5=t(CpGDF1_5)
# CpGDF1_5=data.frame(NbrIndMin=as.numeric(gsub("uniteCov_", "", row.names(CpGDF1_5))),
#                     NbrCpG=CpGDF1_5[,1])
# 
# CpGDF6_12=data.frame(lapply(mylist_uniteCov, function(x) length(x$coverage1)))
# CpGDF6_12=t(CpGDF6_12)
# CpGDF6_12=data.frame(NbrIndMin=as.numeric(gsub("uniteCov_", "", row.names(CpGDF6_12))),
#                      NbrCpG=CpGDF6_12[,1])
# 
# CpGDF=rbind(CpGDF1_5, CpGDF6_12)

# print("Save object for plotting")
# save(CpGDF, file="/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/Data/05MethylKit/plots/CpGDF.RData")

# print("Save plot")
# plotCpGshared <- ggplot(CpGDF, aes(x=NbrIndMin, y=NbrCpG))+
#     geom_smooth(se = F, col = "red")+
#     geom_smooth(method = "lm", se = F, col = "black") +
#     geom_point() +
#     scale_x_continuous("Number of individual fish per treatment group sharing the same methylated CpG sites",
#                        labels = as.character(CpGDF$NbrIndMin), breaks = CpGDF$NbrIndMin)+
#     scale_y_continuous("Number of shared methylated CpG sites") +
#     theme_bw() +
#     geom_hline(yintercept=CpGALL)
# 
# plotCpGshared
# 
# pdf(file="/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/Data/05MethylKit/plots/plotCpGshared.pdf")
# plotCpGshared
# dev.off()
