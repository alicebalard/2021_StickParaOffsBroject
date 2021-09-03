## Methylation analyses - Part 3 
library(plyr) # for join (keep row order)
library(dendextend) # The package dendextend contains many functions for changing the appearance of a dendrogram and for comparing dendrograms.

## load custom functions
source("customRfunctions.R")

## Load previously united data (all 6 treatments)
## uniteCov10: CpG covered in at least 10 individuals (has NAs)
#get(load("/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/Data/05MethylKit/uniteCov10.RData"))

## uniteCovALL: CpG covered in ALL individuals (has no NAs, useful for exploratory clustering analyses)
# get(load("/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/Data/05MethylKit/uniteCovALL.RData"))

## Load samples metadata
metadata144 <- readxl::read_xlsx("/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/GIT_StickParaOffsBroject/data/Kostas_G2_info.xlsx") 
metadata143 <- metadata144[!metadata144$ID %in% "S12",] # S12 was previously removed, utterly bad data
as.factor(metadata144$Family)# 5 families

## Make cluster with colored bars
makePrettyMethCluster(uniteCovALL_mem, metadata143)
## -> cluster by family (as expected, the methylation profile is tightly linked with genes)

## PCA analysis on our samples: plot a scree plot for importance of components
PCASamples(uniteCovALL_mem, screeplot=TRUE)

## Plot PCA for axes 1 & 2
PCASamples(uniteCovALL_mem)
title(sub="colored by family")
## OBVIOUS outliers: S142 & S118

metadata[metadata$ID %in% "S118",]
metadata[metadata$ID %in% "S142",]

################################################################
## creates a new methylBase object by removing the 2 outliers ##
################################################################
uniteCovALL_outrm=reorganize(
  uniteCovALL_mem,
  sample.ids=uniteCovALL_mem@sample.ids[!uniteCovALL_mem@sample.ids %in% c("S118", "S142")],
  treatment=uniteCovALL_mem@treatment[!uniteCovALL_mem@sample.ids %in% c("S118", "S142")])

# update metadata accordingly
metadata141 <- metadata143[!metadata143$ID %in% c("S118", "S142"),]

## Comparison with/without outliers: check that does not change the clustering
makePrettyMethCluster(uniteCovALL_mem, metadata143)
makePrettyMethCluster(uniteCovALL_outrm, metadata141)

## PCA again on several responses, after removal of outliers
PCASamples(uniteCovALL_outrm)
title(sub="colored by family")

## group by SEX
uniteCovALL_outrm_SEX = uniteCovALL_outrm
uniteCovALL_outrm_SEX@treatment = as.numeric(as.factor(metadata141$Sex))
PCASamples(uniteCovALL_outrm_SEX)
title(sub="colored by sex")

## group by treatment
uniteCovALL_outrm_trtG1G2 = uniteCovALL_outrm
uniteCovALL_outrm_trtG1G2@treatment = as.numeric(as.factor(metadata141$trtG1G2))
PCASamples(uniteCovALL_outrm_trtG1G2)
title(sub="colored by trtG1G2")

## group by inf or not
uniteCovALL_outrm_infornot = uniteCovALL_outrm
uniteCovALL_outrm_infornot@treatment = as.numeric(as.factor(metadata141$outcome))
PCASamples(uniteCovALL_outrm_infornot)
title(sub="colored by infected or control")

##############################################################################
## creates a new methylBase object by removing Fam12 (N=4, only in parents) ##
##############################################################################
uniteCovALL_nofam12=reorganize(
  uniteCovALL_outrm,
  sample.ids=uniteCovALL_outrm@sample.ids[
    !uniteCovALL_outrm@sample.ids %in% metadata141$ID[metadata141$Family %in% "Fam12"]],
  treatment=uniteCovALL_outrm@treatment[
    !uniteCovALL_outrm@sample.ids %in% metadata141$ID[metadata141$Family %in% "Fam12"]])

# update metadata accordingly
metadataNoFam12 <- metadata141[!metadata141$Family %in% "Fam12",]

## Comparison with/without outliers: check that does not change the clustering
pdf(file = "~/figTemp/clusterALLCpG.pdf", width = 20, height = 5)
makePrettyMethCluster(uniteCovALL_outrm, metadata141)
dev.off()

pdf(file = "~/figTemp/clusterALLCpG_minusFam12.pdf", width = 20, height = 5)
makePrettyMethCluster(uniteCovALL_nofam12, metadataNoFam12)
dev.off()
