# Each script sources the previous script of the pipeline if needed
source("R05_GlobalMethylationAndFitness.R")

message("R06 starting...\n")
## Produces the subplots for figure 2
### dataOut/fig/Fig2A_clusterALLCpG_offspings.pdf
### dataOut/fig/Fig2B_NMDSplot_allG2.pdf

# Methylation profile, CpG present in all fish
## Dendogram of methylations
print("CpG sites covered in all G1 & G2 fish:")
print(nrow(uniteCovALL_woSexAndUnknowChr))
print("number of fish in total:")
print(length(uniteCovALL_woSexAndUnknowChr@sample.ids))

makePrettyMethCluster(uniteCovALL_woSexAndUnknowChr, fullMetadata,
                      my.cols.trt=c("#333333ff","#ff0000ff",colOffs),
                      my.cols.fam = c(1:4), nbrk = 8)

### Offspring:
print("CpG sites covered in all G2 fish:")
print(nrow(uniteCovALL_G2_woSexAndUnknowChr))
print("number of G2 fish:")
print(length(uniteCovALL_G2_woSexAndUnknowChr@sample.ids))

# Save
pdf(file = "../../dataOut/fig/Fig2A_clusterALLCpG_offspings.pdf", width = 10, height = 4)
makePrettyMethCluster(uniteCovALL_G2_woSexAndUnknowChr, fullMetadata_OFFS,
                      my.cols.fam = c(1:4), nbrk = 8)
dev.off()

## Adonis tests: impact of different variables on methylation pattern
# Permutational multivariate analysis of variance (PERMANOVA) is a
# non-parametric multivariate statistical test. It is used to compare
# groups of objects and test the null hypothesis that the centroids and
# dispersion of the groups as defined by measure space are equivalent for
# all groups.

# make distance matrix with B-C distances
data.dist = makeDatadistFUN(uniteCovALL_G2_woSexAndUnknowChr)

## Adonis test: importance of each predictor
adonis2(data.dist ~ PAT * outcome * Sex * brotherPairID, data = fullMetadata_OFFS, 
        by = "terms")

# Results: family of the father (brotherPairID) explains more than 14% of
# the variance in methylation.
#                                 Df SumOfSqs    R2      F    Pr(>F)    
# PAT                             1 0.002782 0.01470 1.8292  0.001 ***
# outcome                         1 0.001909 0.01009 1.2552  0.027 *  
# Sex                             1 0.002418 0.01277 1.5895  0.001 ***
# brotherPairID                   7 0.028018 0.14805 2.6316  0.001 ***
# PAT:brotherPairID               7 0.014344 0.07580 1.3473  0.001 ***

# To focus on G1 and G2 treatments, we define the permutation structure
# considering brother pairs (N = 8), and use a PERMANOVA to test the
# hypothesis that paternal treatment, offspring treatment and their
# interactions significantly influencing global methylation.
perm <- how(nperm = 1000) # 1000 permutations
setBlocks(perm) <- with(fullMetadata_OFFS, brotherPairID) # define the permutation structure considering brotherPairID and sex
print(adonis2(data.dist ~ PAT * outcome * Sex, data = fullMetadata_OFFS, 
              permutations = perm, by = "terms"))

# Results:
# -   1.5% of the variation explained by PAT (R2=0.01470, p \< 0.001)
# -   1% of the variation explained by outcome (R2=0.01009, p = 0.042)
# -   1.3% of the variation explained by Sex (R2=0.01277, p = 0.011)

## NMDS
#### RUN Goodness of fit
#myGOF.NMDS.FUN(dataset = uniteCovALL_G2_woSexAndUnknowChr)

# Goodness of fit for NMDS suggested the presence of six dimensions with a
# stress value \> 0.1 and 2 with \> 0.2

run = FALSE
if (run==TRUE){
  ## to find the seed that allows convergence:
  # sapply(3:10, function(x) myNMDS(dataset = uniteCovALL_G2_woSexAndUnknowChr, metadata = fullMetadata_OFFS, myseed = x))
  NMDSanalysis <- myNMDSFUN(dataset = uniteCovALL_G2_woSexAndUnknowChr, metadata = fullMetadata_OFFS, myseed = 4)
  pdf(file = "../../dataOut/fig/Fig2B_NMDSplot_allG2.pdf", width = 9, height = 11)
  NMDSanalysis$NMDSplot
  dev.off()
}

## The methylation pattern is more affected by direct treatment when the father was infected (Adonis test WITHIN both parental trt)

### 1. Parents NOT infected
AdonisWithinG1trtFUN(trtgp = c(2,3))

# Results:
# -   brother pair explains 22% of the observed variance p \< 0.01
# -   sex explains 2.3% of the observed variance p \< 0.001
# -   Offspring trt explains 1.6% of the observed variance (non significant)

run = FALSE
if (run==TRUE){
  NMDSanalysis_G1control <- myNMDSFUN(dataset = uniteCovALL_G2_woSexAndUnknowChr,
                                      metadata = fullMetadata_OFFS, myseed = 25,
                                      byParentTrt=TRUE,
                                      trtgp = c(2,3))
  #png(filename = "../../dataOut/NMDSplot_G1fromControlG2.png", width = 900, height = 900)
  NMDSanalysis_G1control$NMDSplot
  #dev.off() 
}

### 2. Parents infected
AdonisWithinG1trtFUN(trtgp = c(5,6))

# Results:
# -   brother pair explains 23% of the observed variance p = 0.01
# -   sex explains 2.1% of the observed variance (non significant)
# -   Offspring trt explains 2.3% of the observed variance p \< 0.01

run = FALSE
if (run==TRUE){
  
  NMDSanalysis_G1infected <- myNMDSFUN(dataset = uniteCovALL_G2_woSexAndUnknowChr,
                                       metadata = fullMetadata_OFFS, myseed = 25,
                                       byParentTrt=TRUE,
                                       trtgp = c(5,6))
  #png(filename = "../../dataOut/NMDSplot_G1fromInfectedG2.png", width = 900, height = 900)
  NMDSanalysis_G1infected$NMDSplot
  #dev.off()
}

message("R06 done\n")
