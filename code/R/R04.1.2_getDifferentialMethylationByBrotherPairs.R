## Differential methylation analyses BY BROTHER PAIR
## A. Balard
## February 2022

#################### Data load & preparation ####################
source("librariesLoading.R")
# load custom functions
source("customRfunctions.R")
## Load samples metadata
source("R02.1_loadMetadata.R")
## define in which machine we're working (apocrita or mythinkpad)
##machine="apocrita"
machine="mythinkpad"
## Load methylation data
loadALL = FALSE # only load CpG shared by half fish per trt group + ALL
source("R02.2_loadMethyldata.R")

######### 
## Calculate parental DMS (parDMS) for each brother pair

## CpG covered in half trt group, parental DMS:
nrow(uniteCov6_G1_woSexAndUnknowChrOVERLAP) # methylBase object with 1001880 rows

## Calculate DMS/DMR
getDiffMethSimple <- function(myuniteCov, myMetadata){
  myDiffMeth=calculateDiffMeth(myuniteCov, mc.cores = 3)#10 on Apocrita
  ## We select the bases that have q-value<0.01 and percent methylation difference larger than 15%.
  ## NB: arg type="hyper" or type="hypo" gives hyper-methylated or hypo-methylated regions/bases.
  myDMS_15pc = getMethylDiff(myDiffMeth, difference=15, qvalue=0.01)
  return(myDMS_15pc)
}

#############
# Select only G1 fish from BP which have offspring
metadata_G1 <- fullMetadata_PAR[
  fullMetadata_PAR$brotherPairID %in% fullMetadata_OFFS$brotherPairID]

#### We will apply thge following function to all BP:
vecBP <- unique(metadata_G2$brotherPairID)
# some parental BP didn't have offspring with sequencing
vecBP <- vecBP[vecBP %in% fullMetadata_OFFS$brotherPairID]

getReacNorm <- function(BP){
  metadata1BP_G1 <- metadata_G1[metadata_G1$brotherPairID %in% BP,]

  # Select methylation object for 1 BP from methylobject with CpG covered in ALL fish
  uniteCov1BP_G1 <- reorganize(methylObj = uniteCovALL_woSexAndUnknowChr,
                               treatment = metadata1BP_G1$trtG1G2_NUM,
                               sample.ids = metadata1BP_G1$ID)
  
  DMS15pc_1BP <- getDiffMethSimple(uniteCov1BP_G1, metadata1BP_G1)
  ## NOTE: performing 'fast.fisher' instead of 'F' for two groups testing.
  
  #############
  ## Get methylation in offspring of these 2 brothers at the parDMS positions:
  
  # metadata
  metadata1BP_G2 <- fullMetadata_OFFS[fullMetadata_OFFS$brotherPairID %in% BP,]
  
  # methyldata
  uniteCov1BP_G2 <- reorganize(methylObj = uniteCovALL_woSexAndUnknowChr,
                               treatment = metadata1BP_G2$trtG1G2_NUM,
                               sample.ids = metadata1BP_G2$ID)
  
  # keep only CpG corresponding to parDMS in this BP, split in hypo/hyper
  parDMS1BP_G2_hypo <- methylKit::select(uniteCov1BP_G2, 
                                         which(paste(uniteCov1BP_G2$chr, uniteCov1BP_G2$start) %in% 
                                                 paste(DMS15pc_1BP[DMS15pc_1BP$meth.diff<0, ]$chr, 
                                                       DMS15pc_1BP[DMS15pc_1BP$meth.diff<0, ]$start)))
  parDMS1BP_G2_hyper <- methylKit::select(uniteCov1BP_G2, 
                                          which(paste(uniteCov1BP_G2$chr, uniteCov1BP_G2$start) %in% 
                                                  paste(DMS15pc_1BP[DMS15pc_1BP$meth.diff>0, ]$chr, 
                                                        DMS15pc_1BP[DMS15pc_1BP$meth.diff>0, ]$start)))
  
  # Get average methylation per individual at these CpGs
  aveMeth_hypo <- apply(percMethylation(parDMS1BP_G2_hypo), 2, mean)
  aveMeth_hyper <- apply(percMethylation(parDMS1BP_G2_hyper), 2, mean)
  
  # make datasets for statistical test
  dfhypoG2parDMS <- data.frame(SampleID=names(aveMeth_hypo), aveMeth=aveMeth_hypo)
  dfhyperG2parDMS <- data.frame(SampleID=names(aveMeth_hyper), aveMeth=aveMeth_hyper)
  
  # add metadata
  dfhypoG2parDMS <- merge(dfhypoG2parDMS, fullMetadata_OFFS)
  dfhyperG2parDMS <- merge(dfhyperG2parDMS, fullMetadata_OFFS)
  
  # and test
  print(paste0("HYPO parDMS for brother pair ", BP))
  modfull <- lmer(aveMeth ~ patTrt * outcome + (1|Sex), data = dfhypoG2parDMS)
  modnointer <- lmer(aveMeth ~ patTrt + outcome + (1|Sex), data = dfhypoG2parDMS)
  modnoG1 <- lmer(aveMeth ~ outcome + (1|Sex), data = dfhypoG2parDMS)
  modnoG2 <- lmer(aveMeth ~ patTrt + (1|Sex), data = dfhypoG2parDMS)
  # test significance of factors
  print("Interactions:")
  print(lrtest(modfull, modnointer))
  print("Parental effect:")
  print(lrtest(modfull, modnoG1))
  print("Offspring effect:")
  print(lrtest(modfull, modnoG2))
  
  print(paste0("HYPER parDMS for brother pair ", BP))
  modfull <- lmer(aveMeth ~ patTrt * outcome + (1|Sex), data = dfhyperG2parDMS)
  modnointer <- lmer(aveMeth ~ patTrt + outcome + (1|Sex), data = dfhyperG2parDMS)
  modnoG1 <- lmer(aveMeth ~ outcome + (1|Sex), data = dfhyperG2parDMS)
  modnoG2 <- lmer(aveMeth ~ patTrt + (1|Sex), data = dfhyperG2parDMS)
  # test significance of factors
  print("Interactions:")
  print(lrtest(modfull, modnointer))
  print("Parental effect:")
  print(lrtest(modfull, modnoG1))
  print("Offspring effect:")
  print(lrtest(modfull, modnoG2))
  
  return(list(dfhypoG2parDMS=dfhypoG2parDMS, dfhyperG2parDMS=dfhyperG2parDMS))
}

results <- sapply(vecBP, getReacNorm) 

results


results[1]


# plot
plist <- list() # empty plot list
for (i in 1:length(results[1:8])){
  plist[[i]] <- ggplot(data = data.frame(results[i]), aes(x=outcome, y=aveMeth, fill=trtG1G2))+
    geom_boxplot()+
    geom_point()+
    scale_fill_manual(values = colOffs) +
    facet_grid(.~patTrt)+
    theme(legend.position = "null", axis.title.x = element_blank()) + 
    ylab("Mean methylation ratio at parental DMS")+
    # scale_y_continuous(limits=c(45,75))
    scale_y_continuous(limits=c(35,75))
} 
do.call(grid.arrange, c(plist, ncol = 4))

## Hyper methylation
plist <- list() # empty plot list
for (i in c(1:8)){
  plist[[i]] <-   ggplot(data = data.frame(results[i+8]), aes(x=outcome, y=aveMeth, fill=trtG1G2))+
    geom_boxplot()+
    geom_point()+
    scale_fill_manual(values = colOffs) +
    facet_grid(.~patTrt)+
    theme(legend.position = "null", axis.title.x = element_blank()) + 
    ylab("Mean methylation ratio at parental DMS")+
    # scale_y_continuous(limits=c(45,75))
    scale_y_continuous(limits=c(35,75))
} 

do.call(grid.arrange, c(plist, ncol = 4))

