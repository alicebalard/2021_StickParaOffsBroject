## MethylKit analysis
## A. Balard
## 25th of August 2021

## R/4.0.2 to run
library(methylKit)
library(readxl)

## Sources:
## https://www.bioconductor.org/packages/devel/bioc/vignettes/methylKit/inst/doc/methylKit.html

# https://nbis-workshop-epigenomics.readthedocs.io/en/latest/content/tutorials/methylationSeq/Seq_Tutorial.html#load-datasets

# /data/archive/archive-SBCS-EizaguirreLab/RRBS/StickPara_Broject_archive/08CompGenomes_mCextr/04RAnalyses_Methylome/methylbam_peichel...

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
## In the following part (R01.2), we decided to remove several samples (N=7)
## S12 (bad quality), S118 & S142 (very weird methylation profiles),
## and Fam12 (N=4, only present in parental group)

IDtoRm= c("S12", "S118", "S142", metadata$ID[metadata$Family %in% "Fam12"])

## create a new methylRawList object
print("Remove 7 samples")
myobj=reorganize(
    myobj,
    sample.ids=metadata$ID[!metadata$ID %in% IDtoRm],
    treatment=metadata$trtG1G2_NUM[!metadata$ID %in% IDtoRm])

## Filtering based on coverage:
# It might be useful to filter samples based on coverage. Particularly, if our samples are suffering from PCR bias it would be useful to discard bases with very high read coverage. Furthermore, we would also like to discard bases that have low read coverage, a high enough read coverage will increase the power of the statistical tests. The code below filters a methylRawList and discards bases that have coverage below 10X and also discards the bases that have more than 99.9th percentile of coverage in each sample.
print("Filter")
filtered.myobj=filterByCoverage(myobj, lo.count=10, lo.perc=NULL,
                                hi.count=NULL, hi.perc=99.9)

## normalise the coverage
print("Normalise")
normFil.myobj=normalizeCoverage(filtered.myobj)

## MERGING SAMPLES: In order to do further analysis, we will need to get the bases covered in all samples. The following function will merge all samples to one object for base-pair locations that are covered in all samples. The unite() function will return a methylBase object which will be our main object for all comparative analysis. The methylBase object contains methylation information for regions/bases that are covered in all samples.

# table(metadata$trtG1G2[!metadata$ID %in% IDtoRm])
#   Control  E_control  E_exposed    Exposed NE_control NE_exposed 
#        12         29         29         12         28         27 

## we kept for downstream analyses all CpG sites present in at least 1 to 12 individuals per group, or in all individuals:
print("Unite and store in a list")
mylist_uniteCov=list()
# or (i in 1:12L){ 
#    uniteCov=unite(normFil.myobj, min.per.group=i, mc.cores=8)# try with 8 cores
#    uniteCov=as(uniteCov,"methylBase")
#    name=paste0("uniteCov_", as.character(i))
#    mylist_uniteCov[[name]]=uniteCov
# 

## Add CpG present in ALL individuals
uniteCov=unite(normFil.myobj, mc.cores=8)
uniteCov=as(uniteCov,"methylBase")
mylist_uniteCov[["uniteCov_ALL"]]=uniteCov

# Idea: plot number of retained CpG site by nbr of individuals sharing these CpG sites. Preparing file for that:
print("Make DF")
CpGDF_ALL=data.frame(NbrIndMin=names(mylist_uniteCov),
                 NbrCpG=lapply(mylist_uniteCov, function(x) length(x$coverage1)))

print("Save object for plotting")
save(CpGDF_ALL, file="/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/Data/05MethylKit/CpGDF_ALL.RData")
