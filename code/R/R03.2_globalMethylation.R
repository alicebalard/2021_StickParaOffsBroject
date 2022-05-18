## Global methylation analyses
## A. Balard
## November 2021

machine="mythinkpad" # define the machine we work on
loadALL = TRUE # load all uniteCov objects
source("R02.3_DATALOAD.R")

#############################################################
### PART 1: Methylation profiles, CpG present in all fish ###
#############################################################

## Dendogram of methylations
## All samples:
pdf("Rfigures/clusterALLCpG.pdf", width = 17, height = 8)
makePrettyMethCluster(uniteCovALL_woSexAndUnknowChr, fullMetadata,
                      my.cols.trt=c("#333333ff","#ff0000ff","#ffe680ff","#ff6600ff","#aaccffff","#aa00d4ff"),
                      my.cols.fam = c(1:4), nbrk = 8)
dev.off()

## offspring: (add TRUE to script loadMethylData)
pdf("Rfigures/clusterALLCpG_offspings.pdf", width = 17, height = 8)
makePrettyMethCluster(uniteCovALL_G2_woSexAndUnknowChr, fullMetadata_OFFS,
                      my.cols.trt=c("#ffe680ff","#ff6600ff", "#aaccffff", "#aa00d4ff"),
                      my.cols.fam = c(1:4), nbrk = 8)
dev.off()

###################
## Let's run Adonis for offspring (considering CpG shared by ALL)
# Permutational multivariate analysis of variance (PERMANOVA) is a non-parametric 
# multivariate statistical test. It is used to compare groups of objects and test
# the null hypothesis that the centroids and dispersion of the groups as defined by
# measure space are equivalent for all groups.
## NB: we define the permutation structure considering brother pairs (N = 8)
table(fullMetadata_OFFS$brotherPairID)

# make distance matrix with B-C distances
data.dist = makeDatadistFUN(uniteCovALL_G2_woSexAndUnknowChr)

## Adonis test: importance of each predictor
adonis2(data.dist ~ PAT * outcome * Sex * brotherPairID, data = fullMetadata_OFFS)
#                                 Df SumOfSqs      R2      F Pr(>F)    
# PAT                             1 0.002782 0.01470 1.8292  0.001 ***
# outcome                         1 0.001909 0.01009 1.2552  0.024 *  
# Sex                             1 0.002418 0.01277 1.5895  0.001 ***
# brotherPairID                   7 0.028018 0.14805 2.6316  0.001 ***
# PAT:brotherPairID               7 0.014344 0.07580 1.3473  0.001 ***

## paternal treatment: explains 1.5% of the variation
## family of the father: 15%
## interaction (individual effect): 7.6%
## sex: 1.3%
## treatment of the offspring: 1%

# We use a PERMANOVA to test the hypothesis that paternal treatment, 
# offspring treatment, sex and their interactions significantly influencing global methylation
perm <- how(nperm = 1000) # 1000 permutations
setBlocks(perm) <- with(fullMetadata_OFFS, brotherPairID) # define the permutation structure considering brotherPairID
## Full model
print(adonis2(data.dist ~ PAT * outcome * Sex, data = fullMetadata_OFFS, permutations = perm))
## remove the non significant interactions
print(adonis2(data.dist ~ PAT + outcome + Sex, data = fullMetadata_OFFS, permutations = perm))
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Blocks:  with(metadata, brotherPairID) 
# Permutation: free
# Number of permutations: 1000
# 
# adonis2(formula = data.dist ~ PAT + outcome + Sex, data = metadata, permutations = perm)
#           Df SumOfSqs      R2      F   Pr(>F)    
# PAT        1 0.002782 0.01470 1.6344 0.000999 *** --> 1.5% of the variation explained by PAT
# outcome    1 0.001909 0.01009 1.1216 0.041958 *  --> 1% of the variation explained by outcome
# Sex        1 0.002418 0.01277 1.4202 0.006993 ** --> 1.3% of the variation explained by Sex
# Residual 107 0.182139 0.96244                    

## Using brotherPairID : G1trt as block
perm <- how(nperm = 1000) # 1000 permutations

dat = fullMetadata_OFFS
dat$brotherPairID_PAT = paste(dat$brotherPairID, dat$PAT)
setBlocks(perm) <- with(dat, brotherPairID_PAT) # define the permutation structure

## with the offspring treatments (result of PAT + outcome)
adonis2(data.dist ~ outcome + Sex, data = dat, permutations = perm)
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Blocks:  with(dat, brotherPairID_PAT) 
# Permutation: free
# Number of permutations: 1000
# 
# adonis2(formula = data.dist ~ outcome + Sex, data = dat, permutations = perm)
# Df SumOfSqs      R2      F   Pr(>F)    
# outcome    1 0.001901 0.01005 1.1113 0.025974 *  <-- within c
# Sex        1 0.002568 0.01357 1.5011 0.000999 ***
# Residual 108 0.184778 0.97638                    
# Total    110 0.189247 1.00000                    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


### Pairwise comparison between treatments:
# myPairAdonisFUN <- function(dataset, metadata){
#   # make distance matrix with B-C distances
#   data.dist = makeDatadistFUN(dataset)
#   pairwise.adonis(x = data.dist, factors = metadata$trtG1G2)
# }
# myPairAdonisFUN(dataset = uniteCovALL_G2_woSexAndUnknowChr, metadata = fullMetadata_OFFS)
## interactions were not significant so we can't expect differences between trtG1G2

########## NMDS
#### RUN Goodness of fit
# myGOF.NMDS.FUN(dataset = uniteCovALL_G2_woSexAndUnknowChr) # Goodness of fit for NMDS 
# suggested the presence of six dimensions with a stress value <0.1 and 2 with > 0.2

## to find the seed that allows convergence: 
# sapply(3:10, function(x) myNMDS(dataset = uniteCovALL_G2_woSexAndUnknowChr, metadata = fullMetadata_OFFS, myseed = x))
NMDSanalysis <- myNMDSFUN(dataset = uniteCovALL_G2_woSexAndUnknowChr, metadata = fullMetadata_OFFS, myseed = 4)
png(filename = "Rfigures/NMDSplot_allG2.png", width = 900, height = 1100)
NMDSanalysis$NMDSplot
dev.off()

#############################################################################################################
## Adonis test WITHIN both parental trt: are G2 from infected G1 more homogeneous than G2 from control G1? ##
#############################################################################################################
table(fullMetadata$trtG1G2, fullMetadata$trtG1G2_NUM)

## 1. Parents are NOT infected. Is there a difference in global methylation between infected/control trt?
AdonisWithinG1trtFUN(trtgp = c(2,3))
#               Df SumOfSqs      R2      F   Pr(>F)   
# outcome        1 0.001475 0.01622 1.0073 0.509491   
# Sex            1 0.002094 0.02302 1.4301 0.002997 **
# brotherPairID  7 0.020026 0.22018 1.9537 0.003996 **
## Clutch explains 22% of the observed variance
## Offspring trt explains 1.6% of the observed variance (non significant)

NMDSanalysis_G1control <- myNMDSFUN(dataset = uniteCovALL_G2_woSexAndUnknowChr, 
                                    metadata = fullMetadata_OFFS, myseed = 25,
                                    byParentTrt=TRUE,
                                    trtgp = c(2,3))
png(filename = "Rfigures/NMDSplot_G1fromControlG2.png", width = 900, height = 900)
NMDSanalysis_G1control$NMDSplot
dev.off()

## 2. Parents are infected. Is there a difference in global methylation between infected/control trt?
AdonisWithinG1trtFUN(trtgp = c(5,6))
#               Df SumOfSqs      R2      F   Pr(>F)   
# outcome        1 0.002160 0.02262 1.4047 0.003996 **
# Sex            1 0.002053 0.02150 1.3349 0.174825   
# brotherPairID  7 0.022076 0.23117 2.0507 0.018981 * 
## Clutch explains 23% of the observed variance
## Offspring trt explains 2.3% of the observed variance


# testConvergence <-sapply(1:20, function(x){
#   run = myNMDSFUN(dataset = uniteCovALL_G2_woSexAndUnknowChr,
#                   metadata = fullMetadata_OFFS, myseed = x,
#                   byParentTrt=TRUE,
#                   trtgp = c(5,6))
#   return(run$NMDS$converged)})# 
# which(testConvergence %in% TRUE)

NMDSanalysis_G1infected <- myNMDSFUN(dataset = uniteCovALL_G2_woSexAndUnknowChr, 
                                     metadata = fullMetadata_OFFS, myseed = 10,
                                     byParentTrt=TRUE,
                                     trtgp = c(5,6))
png(filename = "Rfigures/NMDSplot_G1fromInfectedG2.png", width = 900, height = 900)
NMDSanalysis_G1infected$NMDSplot
dev.off()

##################################################################################################
## Is the BC methylation distance between clutches bigger between BROTHER PAIR or PATERNAL TRT? ##
##################################################################################################

percMethMatG2 = makePercentMetMat(uniteCovALL_G2_woSexAndUnknowChr)
## Percentage methylation for:
nrow(percMethMatG2) # 111 samples
ncol(percMethMatG2) # 78246 positions

percMethMatG2_controlP <- percMethMatG2[
  rownames(percMethMatG2) %in% fullMetadata_OFFS$SampleID[fullMetadata_OFFS$patTrt %in% "controlP"],]
nrow(percMethMatG2_controlP)# 55 samples

percMethMatG2_infectedP <- percMethMatG2[
  rownames(percMethMatG2) %in% fullMetadata_OFFS$SampleID[fullMetadata_OFFS$patTrt %in% "infectedP"],]
nrow(percMethMatG2_infectedP)# 56 samples

## Multiple Response Permutation Procedure (MRPP) provides a test of whether there 
## is a significant difference between two or more groups of sampling units:
G2.mrpp_controlP <- with(fullMetadata_OFFS[fullMetadata_OFFS$patTrt %in% "controlP", ],
                mrpp(percMethMatG2_controlP, Family, distance = "bray"))
G2.mrpp_controlP

G2.mrpp_infectedP <- with(fullMetadata_OFFS[fullMetadata_OFFS$patTrt %in% "infectedP", ],
                         mrpp(percMethMatG2_infectedP, Family, distance = "bray"))
G2.mrpp_infectedP

## Multiple Response Permutation Procedure (MRPP) provides a test of whether there 
## is a significant difference between two or more groups of sampling units:
G2.mrpp <- with(fullMetadata_OFFS, mrpp(percMethMatG2, clutch.ID, distance = "bray"))
G2.mrpp

G2.mrpp_brotherPairID <- with(fullMetadata_OFFS, mrpp(percMethMatG2, brotherPairID, distance = "bray"))
G2.mrpp_brotherPairID

G2.mrpp_patTrt <- with(fullMetadata_OFFS, mrpp(percMethMatG2, patTrt, distance = "bray"))
G2.mrpp_patTrt

## meandist between each of the 16 clutches
dat = fullMetadata_OFFS
dat$trtG1G2 <- gsub("exposed", "E", dat$trtG1G2)
dat$trtG1G2 <- gsub("control", "C", dat$trtG1G2)
dat$trtG1G2 <- gsub("NE", "C", dat$trtG1G2)
dat$clutch.ID <- paste(dat$trtG1G2, sapply(strsplit(as.character(dat$clutch.ID), "_"), `[`, 2))

G2.md <- with(dat, meandist(vegdist(percMethMatG2), clutch.ID, distance = "bray"))
G2.md
summary(G2.md)
heatmap(G2.md)
fullMetadata_OFFS$clutch.ID


# fullMetadata_OFFS$clutch.ID
# fullMetadata_OFFS$brotherPairID
# fullMetadata_OFFS$patTrt
