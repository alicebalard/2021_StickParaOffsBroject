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

⎄###############################################################
## Remove reads from sex chromosome X ("Gy_chrXIX") and unmapped contigs ("Gy_chrUn")

length(uniteCovALL$chr) ## nbr CpG shared by all 137 samples: 47238

length(uniteCov6$chr) ## nbr CpG shared by at least 6 animals: 1319604

## nbr CpGper chrom shared by all 137 samples:
table(uniteCovALL$chr)
## Gy_chrI    Gy_chrII   Gy_chrIII    Gy_chrIV    Gy_chrIX    Gy_chrUn     Gy_chrV    Gy_chrVI   Gy_chrVII  Gy_chrVIII     Gy_chrX 
# 3147        2305        2088        3263        2171        2581        1983        1831        3138        1929        1892 
# Gy_chrXI   Gy_chrXII  Gy_chrXIII   Gy_chrXIV   Gy_chrXIX    Gy_chrXV   Gy_chrXVI  Gy_chrXVII Gy_chrXVIII    Gy_chrXX   Gy_chrXXI 
# 2323        2149        2085        1897         918        2258        1857        2127        1863        1714        1719 

## nbr CpG on sex chromosome of unmapped: 3499
nrow(uniteCovALL[uniteCovALL$chr %in% c("Gy_chrXIX", "Gy_chrUn"),])

## Keep CpG apart from sex chromosome XIX and unmapped (comprise Y chr)
## CpG OUTSIDE of these chromosomes:43739
uniteCovALL_final=uniteCovALL[!uniteCovALL$chr %in% c("Gy_chrXIX", "Gy_chrUn"),]
nrow(uniteCovALL_final)

## CpG  
uniteCov6_final=uniteCov6[!uniteCov6$chr %in% c("Gy_chrXIX", "Gy_chrUn"),]
nrow(uniteCov6_final) ## 1183060

#######################################
## Saving point for further analyses ##
#######################################
save(uniteCovALL_final, file = "/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/Data/05MethylKit/uniteCovALL_final_woSexAndUnknownChr.RData")
save(uniteCov6_final, file = "/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/Data/05MethylKit/uniteCov6_final_woSexAndUnknownChr.RData")

## Comparison with/without outliers: check that does not change the clustering
pdf(file = "/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/GIT_StickParaOffsBroject/data/fig/clusterALLCpG.pdf", 
    width = 20, height = 7)
makePrettyMethCluster(uniteCovALL_N137_final, metadata)
dev.off()

################# PCA
## PCA analysis on our samples: plot a scree plot for importance of components
p=PCASamples(uniteCovALL_N137_final, screeplot=TRUE, obj.return = T) # first axis very important
s=summary(p)
#create scree plot
library(ggplot2)
qplot(c(1:137), s$importance[2,]) + 
  geom_line() + 
  xlab("Principal Component") + 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1L), name = "Variance Explained")+
  ggtitle("Scree Plot") +
  theme_bw()

## PCA: validate the absence of strong sex influence on the methylation
uniteCovALL_N137_final_SEX = uniteCovALL_N137_final
uniteCovALL_N137_final_SEX@treatment = as.numeric(as.factor(metadata$Sex))
PCASamples(uniteCovALL_N137_final_SEX); title(sub="colored by sex")
PCASamples(uniteCovALL_N137_final_SEX, comp = c(3,4)); title(sub="colored by sex")

## PCA: check the Family influence on the methylation pattern
uniteCovALL_N137_final_FAM = uniteCovALL_N137_final
uniteCovALL_N137_final_FAM@treatment = as.numeric(as.factor(metadata$Family))
PCASamples(uniteCovALL_N137_final_FAM); title(sub="colored by family")
PCASamples(uniteCovALL_N137_final_FAM, comp = c(3,4)); title(sub="colored by family")

## PCA: check the Pattern treatment influence on the methylation pattern
uniteCovALL_N137_final_PAT = uniteCovALL_N137_final
metadata$PAT="Exposed father group"
metadata$PAT[metadata$trtG1G2 %in% c("Control", "NE_control", "NE_exposed")]="Control father group"
uniteCovALL_N137_final_PAT@treatment = as.numeric(as.factor(metadata$PAT))
PCASamples(uniteCovALL_N137_final_PAT); title(sub="colored by paternal exposure group")
PCASamples(uniteCovALL_N137_final_PAT, comp=c(3,4)); title(sub="colored by paternal exposure group")

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
