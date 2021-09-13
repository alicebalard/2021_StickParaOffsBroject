## Methylation analyses - Part 3
library(plyr) # for join (keep row order)
library(dendextend) # The package dendextend contains many functions for changing the appearance of a dendrogram and for comparing dendrograms.
library(ggplot2)

## load custom functions
source("/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/GIT_StickParaOffsBroject/code/R/customRfunctions.R")

## Load previously united data (all 6 treatments)
## uniteCovALL: CpG covered in ALL individuals (has no NAs, useful for exploratory clustering analyses)
get(load("/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/Data/05MethylKit/uniteCovALL_N137.RData"))
## For further analyses: CpG covered in at least 6 individuals per group
get(load("/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/Data/05MethylKit/uniteCov6_N137.RData"))

## Load samples metadata
metadata144 <- readxl::read_xlsx("/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/GIT_StickParaOffsBroject/data/Kostas_G2_info.xlsx")
metadata144$trtG1G2_NUM <- as.numeric(as.factor(metadata144$trtG1G2))
## remove 7 samples (3 outliers + 1 family of 4 only present in parents)
IDtoRm= c("S12", "S118", "S142", metadata144$ID[metadata144$Family %in% "Fam12"])
metadata <- metadata144[!metadata144$ID %in% IDtoRm,] 

#### To be updated with final one and that chunk will be removed
## create a new methylRawList object
uniteCovALL=reorganize(
  uniteCovALL_N137,
  sample.ids=metadata$ID,
  treatment=metadata$trtG1G2_NUM)

uniteCov6=reorganize(
  uniteCov6_N137,
  sample.ids=metadata$ID,
  treatment=metadata$trtG1G2_NUM)

###############################################################
## Remove reads from sex chromosome X ("Gy_chrXIX") and unmapped contigs ("Gy_chrUn")

## nbr CpG=47238
length(uniteCovALL_N137$chr)

## nbr CpGper chrom=
# Gy_chrI    Gy_chrII   Gy_chrIII    Gy_chrIV    Gy_chrIX    Gy_chrUn     Gy_chrV    Gy_chrVI   Gy_chrVII  Gy_chrVIII     Gy_chrX 
# 3147        2305        2088        3263        2171        2581        1983        1831        3138        1929        1892 
# Gy_chrXI   Gy_chrXII  Gy_chrXIII   Gy_chrXIV   Gy_chrXIX    Gy_chrXV   Gy_chrXVI  Gy_chrXVII Gy_chrXVIII    Gy_chrXX   Gy_chrXXI 
# 2323        2149        2085        1897         918        2258        1857        2127        1863        1714        1719 
table(uniteCovALL_N137$chr)

## nbr CpG on sex chromosome of unmapped: 3499
uniteCovALL_N137[uniteCovALL_N137$chr %in% c("Gy_chrXIX", "Gy_chrUn"),]

## CpG OUTSIDE of these chromosomes:43739
uniteCovALL_N137_final=uniteCovALL_N137[!uniteCovALL_N137$chr %in% c("Gy_chrXIX", "Gy_chrUn"),]

## Comparison with/without outliers: check that does not change the clustering
pdf(file = "/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/GIT_StickParaOffsBroject/data/fig/clusterALLCpG.pdf", 
    width = 20, height = 7)
makePrettyMethCluster(uniteCovALL_N137_final, metadata)
dev.off()


###############################################################
## When at least 6 individuals per group share the CpGs:1183060
uniteCov6_N137_final=uniteCov6_N137[!uniteCov6_N137$chr %in% c("Gy_chrXIX", "Gy_chrUn"),]
length(uniteCov6_N137_final$chr)





####### Previous code if needed to repeat for supplements
# ## Load samples metadata
# metadata144 <- readxl::read_xlsx("/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/GIT_StickParaOffsBroject/data/Kostas_G2_info.xlsx") 
# metadata143 <- metadata144[!metadata144$ID %in% "S12",] # S12 was previously removed, utterly bad data
# as.factor(metadata144$Family)# 5 families
# 
# ## Make cluster with colored bars
# makePrettyMethCluster(uniteCovALL_mem, metadata143)
# ## -> cluster by family (as expected, the methylation profile is tightly linked with genes)
# 
# ## PCA analysis on our samples: plot a scree plot for importance of components
# PCASamples(uniteCovALL_mem, screeplot=TRUE)
# 
# ## Plot PCA for axes 1 & 2
# PCASamples(uniteCovALL_mem)
# title(sub="colored by family")
# ## OBVIOUS outliers: S142 & S118
# 
# metadata[metadata$ID %in% "S118",]
# metadata[metadata$ID %in% "S142",]
# 
# ################################################################
# ## creates a new methylBase object by removing the 2 outliers ##
# ################################################################
# uniteCovALL_outrm=reorganize(
#   uniteCovALL_mem,
#   sample.ids=uniteCovALL_mem@sample.ids[!uniteCovALL_mem@sample.ids %in% c("S118", "S142")],
#   treatment=uniteCovALL_mem@treatment[!uniteCovALL_mem@sample.ids %in% c("S118", "S142")])
# 
# # update metadata accordingly
# metadata141 <- metadata143[!metadata143$ID %in% c("S118", "S142"),]
# 
# ## Comparison with/without outliers: check that does not change the clustering
# makePrettyMethCluster(uniteCovALL_mem, metadata143)
# makePrettyMethCluster(uniteCovALL_outrm, metadata141)
# 
# ## PCA again on several responses, after removal of outliers
# PCASamples(uniteCovALL_outrm)
# title(sub="colored by family")
# 
# ## group by SEX
# uniteCovALL_outrm_SEX = uniteCovALL_outrm
# uniteCovALL_outrm_SEX@treatment = as.numeric(as.factor(metadata141$Sex))
# PCASamples(uniteCovALL_outrm_SEX)
# title(sub="colored by sex")
# 
# ## group by treatment
# uniteCovALL_outrm_trtG1G2 = uniteCovALL_outrm
# uniteCovALL_outrm_trtG1G2@treatment = as.numeric(as.factor(metadata141$trtG1G2))
# PCASamples(uniteCovALL_outrm_trtG1G2)
# title(sub="colored by trtG1G2")
# 
# ## group by inf or not
# uniteCovALL_outrm_infornot = uniteCovALL_outrm
# uniteCovALL_outrm_infornot@treatment = as.numeric(as.factor(metadata141$outcome))
# PCASamples(uniteCovALL_outrm_infornot)
# title(sub="colored by infected or control")
# 
# ##############################################################################
# ## creates a new methylBase object by removing Fam12 (N=4, only in parents) ##
# ##############################################################################
# uniteCovALL_nofam12=reorganize(
#   uniteCovALL_outrm,
#   sample.ids=uniteCovALL_outrm@sample.ids[
#     !uniteCovALL_outrm@sample.ids %in% metadata141$ID[metadata141$Family %in% "Fam12"]],
#   treatment=uniteCovALL_outrm@treatment[
#     !uniteCovALL_outrm@sample.ids %in% metadata141$ID[metadata141$Family %in% "Fam12"]])
# 
# # update metadata accordingly
# metadataNoFam12 <- metadata141[!metadata141$Family %in% "Fam12",]
# 
# ## Comparison with/without outliers: check that does not change the clustering
# pdf(file = "~/figTemp/clusterALLCpG.pdf", width = 20, height = 5)
# makePrettyMethCluster(uniteCovALL_outrm, metadata141)
# dev.off()
# 
# pdf(file = "~/figTemp/clusterALLCpG_minusFam12.pdf", width = 20, height = 5)
# makePrettyMethCluster(uniteCovALL_nofam12, metadataNoFam12)
# dev.off()
