## Global methylation analyses
## A. Balard
## November 2021

#################### Data load & preparation ####################
source("librariesLoading.R")
## load custom functions
source("customRfunctions.R")
## Load samples metadata
source("R02.1_loadMetadata.R")
## define in which machine we're working (apocrita or mythinkpad)
#machine="apocrita"
machine="mythinkpad"
## Load methylation data
loadALL = TRUE # load all uniteCov objects
source("R02.2_loadMethyldata.R")

#############################################################
### PART 1: Methylation profiles, CpG present in all fish ###
#############################################################

## Dendogram of methylations
## All samples:
pdf("Rfigures/clusterALLCpG.pdf", width = 16, height = 7)
makePrettyMethCluster(uniteCovALL_woSexAndUnknowChr, fullMetadata,
                      my.cols.trt=c("#333333ff","#ff0000ff","#ffe680ff","#ff6600ff","#aaccffff","#aa00d4ff"),
                      my.cols.fam = c(1:4))
dev.off()

## offspring: (add TRUE to script loadMethylData)
pdf("Rfigures/clusterALLCpG_offspings.pdf", width = 17, height = 7)
makePrettyMethCluster(uniteCovALL_G2_woSexAndUnknowChr, fullMetadata_OFFS,
                      my.cols.trt=c("#ffe680ff","#ff6600ff", "#aaccffff", "#aa00d4ff"),
                      my.cols.fam = c(1:4))
dev.off()

###################
## Let's run Adonis for offspring (considering CpG shared by ALL)
# Permutational multivariate analysis of variance (PERMANOVA) is a non-parametric 
# multivariate statistical test. It is used to compare groups of objects and test
# the null hypothesis that the centroids and dispersion of the groups as defined by
# measure space are equivalent for all groups.
myadonisFUN(dataset = uniteCovALL_G2_woSexAndUnknowChr, metadata = fullMetadata_OFFS)
## Interactions not significant - run without:
# adonis2(formula = data.dist ~ PAT + outcome + Sex, data = metadata, permutations = perm)
#            Df SumOfSqs      R2      F   Pr(>F)    
# PAT        1 0.002782 0.01470 1.6344 0.000999 *** --> 1.5% of the variation explained by PAT
# outcome    1 0.001909 0.01009 1.1216 0.070929 .  --> 1% of the variation explained by outcome
# Sex        1 0.002418 0.01277 1.4202 0.004995 ** --> 1.3% of the variation explained by Sex
# Residual 107 0.182139 0.96244                    
# Total    110 0.189247 1.00000                    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Blocks:  with(metadata, Family) 
# Permutation: free
# Number of permutations: 1000

######################################################################
## Test the effect of treatment and sex in both groups by parental trt
## 1. Parents are NOT infected. Is there a difference in global methylation between infected/control trt?
trtgp = c(2,3)
# make distance matrix with B-C distances
data.dist = makeDatadistFUN(reorganize(methylObj = uniteCovALL_G2_woSexAndUnknowChr,
                                       treatment = fullMetadata_OFFS$trtG1G2_NUM[fullMetadata_OFFS$trtG1G2_NUM %in% trtgp],
                                       sample.ids = fullMetadata_OFFS$ID[fullMetadata_OFFS$trtG1G2_NUM %in% trtgp]))
perm <- how(nperm = 1000) # 1000 permutations
setBlocks(perm) <- with(fullMetadata_OFFS[fullMetadata_OFFS$trtG1G2_NUM %in% trtgp,], Family) # define the permutation structure considering family
## Full model
adonis2(data.dist ~ outcome * Sex, data = fullMetadata_OFFS[fullMetadata_OFFS$trtG1G2_NUM %in% trtgp,], permutations = perm)
## remove the non significant interactions
adonis2(data.dist ~ outcome + Sex, data = fullMetadata_OFFS[fullMetadata_OFFS$trtG1G2_NUM %in% trtgp,], permutations = perm)
## only sex is significant! It explains 2.3% of variation (outcome: 1.6%, non signif)

## 2. Parents are infected. Is there a difference in global methylation between infected/control trt?
# make distance matrix with B-C distances
trtgp = c(5,6)
data.dist = makeDatadistFUN(reorganize(methylObj = uniteCovALL_G2_woSexAndUnknowChr,
                                       treatment = fullMetadata_OFFS$trtG1G2_NUM[fullMetadata_OFFS$trtG1G2_NUM %in% trtgp],
                                       sample.ids = fullMetadata_OFFS$ID[fullMetadata_OFFS$trtG1G2_NUM %in% trtgp]))
perm <- how(nperm = 1000) # 1000 permutations
setBlocks(perm) <- with(fullMetadata_OFFS[fullMetadata_OFFS$trtG1G2_NUM %in% trtgp,], Family) # define the permutation structure considering family
## Full model
adonis2(data.dist ~ outcome * Sex, data = fullMetadata_OFFS[fullMetadata_OFFS$trtG1G2_NUM %in% trtgp,], permutations = perm)
## remove the non significant interactions
adonis2(data.dist ~ outcome + Sex, data = fullMetadata_OFFS[fullMetadata_OFFS$trtG1G2_NUM %in% trtgp,], permutations = perm)
## outcome is significant! It explains 2.3% of variation (Sex: 2.2%, non signif)

########## NMDS
#### RUN Goodness of fit
# myGOF.NMDS.FUN(dataset = uniteCovALL_G2_woSexAndUnknowChr) # Goodness of fit for NMDS 
# suggested the presence of six dimensions with a stress value <0.1 and 2 with > 0.2
NMDSanalysis <- myNMDS(dataset = uniteCovALL_G2_woSexAndUnknowChr, metadata = fullMetadata_OFFS)
NMDSanalysis$NMDSplot
save(NMDSanalysis, file = "../../data/fig/NMDSplots.RData")
