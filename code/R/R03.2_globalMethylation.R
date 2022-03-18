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
myadonisFUN(dataset = uniteCovALL_G2_woSexAndUnknowChr, metadata = fullMetadata_OFFS)
# PAT        1 0.002782 0.01470 1.6344 0.000999 ***
#   outcome    1 0.001909 0.01009 1.1216 0.068931 .  
# Sex        1 0.002418 0.01277 1.4202 0.010989 *  

# trtG1G2    3 0.006418 0.03391 1.2568 0.000999 ***
#   Sex        1 0.002407 0.01272 1.4141 0.008991 ** 

########## NMDS
#### RUN Goodness of fit
# myGOF.NMDS.FUN(dataset = uniteCovALL_G2_woSexAndUnknowChr) # Goodness of fit for NMDS 
# suggested the presence of six dimensions with a stress value <0.1 and 2 with > 0.2
NMDSanalysis <- myNMDS(dataset = uniteCovALL_G2_woSexAndUnknowChr, metadata = fullMetadata_OFFS)
NMDSanalysis$NMDSplot
save(NMDSanalysis, file = "../../data/fig/NMDSplots.RData")
