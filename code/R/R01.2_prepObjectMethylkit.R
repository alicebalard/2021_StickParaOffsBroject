## MethylKit object preparation
## A. Balard
## 25th of August 2021

## NB: change R01.4 loadMethyldata each time this script is amended

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
source("customRfunctions.R")

##### Load prepared dataset (in APOCRITA) #####
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
uniteCovALL= methylKit::unite(normFil.myobj, mc.cores=8)
uniteCovALL=as(uniteCovALL,"methylBase")

## In all PARENTS:
uniteCovALL_G1 = reorganize(
    normFil.myobj,
    sample.ids=metadata$ID[!metadata$ID %in% IDtoRm &
                           metadata$Generat %in% "P"],
    treatment=metadata$trtG1G2_NUM[!metadata$ID %in% IDtoRm &
                                   metadata$Generat %in% "P"])

uniteCovALL_G1 = methylKit::unite(uniteCovALL_G1, mc.cores=8)# try with 8 cores
uniteCovALL_G1 = as(uniteCovALL_G1,"methylBase")

## In all OFFSPRING:
uniteCovALL_G2 = reorganize(
    normFil.myobj,
    sample.ids=metadata$ID[!metadata$ID %in% IDtoRm &
                           metadata$Generat %in% "O"],
    treatment=metadata$trtG1G2_NUM[!metadata$ID %in% IDtoRm &
                                   metadata$Generat %in% "O"])

uniteCovALL_G2 = methylKit::unite(uniteCovALL_G2, mc.cores=8)# try with 8 cores
uniteCovALL_G2 = as(uniteCovALL_G2,"methylBase")

## Keep methylated CpG sites observed in at least 50% individual fish after filtering and normalising: 6 for parents, 14 for offsprings
## PARENTS
uniteCov6_G1 = reorganize(
    normFil.myobj,
    sample.ids=metadata$ID[!metadata$ID %in% IDtoRm &
                           metadata$Generat %in% "P"],
    treatment=metadata$trtG1G2_NUM[!metadata$ID %in% IDtoRm &
                                   metadata$Generat %in% "P"])

uniteCov6_G1 = methylKit::unite(uniteCov6_G1,
                     min.per.group=6L, mc.cores=8)# try with 8 cores
uniteCov6_G1 = as(uniteCov6_G1,"methylBase")

## OFFSPRING
uniteCov14_G2 = reorganize(
    normFil.myobj,
    sample.ids=metadata$ID[!metadata$ID %in% IDtoRm &
                           metadata$Generat %in% "O"],
    treatment=metadata$trtG1G2_NUM[!metadata$ID %in% IDtoRm &
                                   metadata$Generat %in% "O"])

uniteCov14_G2 = methylKit::unite(uniteCov14_G2,
                     min.per.group=14L, mc.cores=8)# try with 8 cores
uniteCov14_G2 = as(uniteCov14_G2,"methylBase")

## And save in RData file for later analyses
save(uniteCovALL, file = "/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/Data/05MethylKit/output/uniteCovALL.RData")
save(uniteCovALL_G1, file = "/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/Data/05MethylKit/output/uniteCovALL_G1.RData")
save(uniteCovALL_G2, file = "/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/Data/05MethylKit/output/uniteCovALL_G2.RData")
save(uniteCov6_G1, file = "/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/Data/05MethylKit/output/uniteCov6_G1.RData")
save(uniteCov14_G2, file = "/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/Data/05MethylKit/output/uniteCov14_G2.RData")

########################################################################################
## Remove reads from sex chromosome X ("Gy_chrXIX") and unmapped contigs ("Gy_chrUn") ##
########################################################################################
print("nbr CpG shared by all 135 samples:")
length(uniteCovALL$chr)
print("nbr CpG shared by all parents:")
length(uniteCovALL_G1$chr)
print("nbr CpG shared by all offsprings:")
length(uniteCovALL_G2$chr)
print("nbr CpG shared by 50% (6) of parents per trt group:")
length(uniteCov6_G1$chr)
print("nbr CpG shared by 50% (14) of offsprings per trt group:")
length(uniteCov14_G2$chr)

print("nbr CpG on sex chromosome of unmapped:")
nrow(uniteCovALL[uniteCovALL$chr %in% c("Gy_chrXIX", "Gy_chrUn"),])

## Keep CpG apart from sex chromosome XIX and unmapped (comprise Y chr)
uniteCovALL_woSexAndUnknowChr=uniteCovALL[
    !uniteCovALL$chr %in% c("Gy_chrXIX", "Gy_chrUn"),]
uniteCovALL_G1_woSexAndUnknowChr=uniteCovALL_G1[
    !uniteCovALL_G1$chr %in% c("Gy_chrXIX", "Gy_chrUn"),]
uniteCovALL_G2_woSexAndUnknowChr=uniteCovALL_G2[
    !uniteCovALL_G2$chr %in% c("Gy_chrXIX", "Gy_chrUn"),]
uniteCov6_G1_woSexAndUnknowChr=uniteCov6_G1[
    !uniteCov6_G1$chr %in% c("Gy_chrXIX", "Gy_chrUn"),]
uniteCov14_G2_woSexAndUnknowChr=uniteCov14_G2[
    !uniteCov14_G2$chr %in% c("Gy_chrXIX", "Gy_chrUn"),]

#######################################
## Saving point for further analyses ##
#######################################
save(uniteCovALL_woSexAndUnknowChr, file = "/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/Data/05MethylKit/output/uniteCovALL_woSexAndUnknownChr.RData")
save(uniteCovALL_G1_woSexAndUnknowChr, file = "/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/Data/05MethylKit/output/uniteCovALL_G1_woSexAndUnknownChr.RData")
save(uniteCovALL_G2_woSexAndUnknowChr, file = "/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/Data/05MethylKit/output/uniteCovALL_G2_woSexAndUnknownChr.RData")
save(uniteCov6_G1_woSexAndUnknowChr, file = "/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/Data/05MethylKit/output/uniteCov6_G1_woSexAndUnknownChr.RData")
save(uniteCov14_G2_woSexAndUnknowChr, file = "/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/Data/05MethylKit/output/uniteCov14_G2_woSexAndUnknownChr.RData")

## Put the correct date in bash:
# for f in uniteCov*; do mv "$f" "$(echo "$f" | sed s/.RData/_10feb22.RData/)"; done

######################################################################################################
## consider only methylKit object (uniteCov) with OVERLAPPING positions in Parents G1 and Offspring G2
overlappingCpG_G1G2df <- findOverlaps(as(uniteCov6_G1_woSexAndUnknowChr,"GRanges"), 
                                      as(uniteCov14_G2_woSexAndUnknowChr,"GRanges"))
overlappingCpG_G1G2df <- data.frame(overlappingCpG_G1G2df)

uniteCov6_G1_woSexAndUnknowChrOVERLAP <- uniteCov6_G1_woSexAndUnknowChr[overlappingCpG_G1G2df$queryHits,]
uniteCov14_G2_woSexAndUnknowChrOVERLAP <- uniteCov14_G2_woSexAndUnknowChr[overlappingCpG_G1G2df$subjectHits,]

## Save in 2 bits
save(uniteCov6_G1_woSexAndUnknowChrOVERLAP,
     file = "/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/Data/05MethylKit/output/uniteCov6_G1_woSexAndUnknowChrOVERLAP_23feb2022.RData")

save(uniteCov14_G2_woSexAndUnknowChrOVERLAP,
     file = "/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/Data/05MethylKit/output/uniteCov14_G2_woSexAndUnknowChrOVERLAP_23feb2022.RData")
    
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
