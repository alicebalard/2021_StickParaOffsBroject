## Methylation analyses - Part 4
library(methylKit)

###################################################################
## Load methylkit raw object before filtering, normalising, uniting:
load(file = "/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/Data/05MethylKit/myMethylKitObj_BeforeFiltering_137.Rdata")

## Load data from 137 samples, sex chromosome XIX and unmapped chr removed
### First, CpG shared by all individuals: uniteCovALL_final
get(load(file = "/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/Data/05MethylKit/uniteCovALL_final_woSexAndUnknownChr.RData"))

## Then, CpG shared by at least 6 individuals: uniteCov6_final
get(load(file = "/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/Data/05MethylKit/uniteCov6_final_woSexAndUnknownChr.RData"))

## Load metadata file 
metadata <- readxl::read_xlsx("/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/GIT_StickParaOffsBroject/data/Kostas_G2_info.xlsx")
metadata$trtG1G2_NUM <- as.numeric(as.factor(metadata$trtG1G2))
## We removed several samples (N=7)
## S12 (bad quality), S118 & S142 (very weird methylation profiles),
## and Fam12 (N=4, only present in parental group)
IDtoRm= c("S12", "S118", "S142", metadata$ID[metadata$Family %in% "Fam12"])
metadata_137 <- metadata[!metadata$ID %in% IDtoRm, ]
####################################################################

## Question: is there a difference in methylation between offsprings from the same treatment group but with parents from one or the other treatment?
## test:



uniteCov6_final@treatment

metadata_137$trtG1G2
metadata_137$trtG1G2_NUM


nrow(uniteCov6_final)

1+1
