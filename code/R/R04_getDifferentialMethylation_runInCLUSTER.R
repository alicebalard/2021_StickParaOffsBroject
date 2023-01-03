## Differential methylation analyses
## A. Balard
## February 2022 (updated Dec22)

# Each script sources the previous script of the pipeline if needed
source("R03_prepObjectMethylkit_runInCLUSTER.R")

###############################################
### Differential Methylation Sites/Regions ####
###############################################

### 5 comparisons to do:
# 0. PARENTS trt-ctrl
# 1. CC-TC = CONTROL fish (parent CvsT)
# 2. CT-TT = TREATMENT fish (parent CvsT)
# 3. CC-CT = fish from CONTROL parents (G2 CvsT)
# 4. TC-TT = fish from TREATMENT parents (G2 CvsT)

## 1&2: paternal effect 
## 3&4: G2 direct treatment effect

################################################################
## Correspondance trtG1G2 and numerical values used by MethylKit
# table(fullMetadata$trtG1G2_NUM, fullMetadata$trtG1G2)
## Control C 1
## Exposed T 4
## NE_control CC 5
## NE_exposed CT 6
## E_control TC 2
## E_exposed TT 3

#####################################################################
## Calculate DMS/DMR accounting for covariates: brotherPairID and sex

# 0. PARENTS trt-ctrl
# 1. CC-TC = CONTROL fish (parent CvsT)
# 2. CT-TT = TREATMENT fish (parent CvsT)
# 3. CC-CT = fish from CONTROL parents (G2 CvsT)
# 4. TC-TT = fish from TREATMENT parents (G2 CvsT)

## Per brother pair:
getDMperBP <- function(BP){
  ## Unite object for one Brother Pair:
  metadataBP_CC_TC = fullMetadata_OFFS[fullMetadata_OFFS$brotherPairID %in% BP &
                                         fullMetadata_OFFS$trtG1G2 %in% c("NE_control", "E_control"), ]
  metadataBP_CT_TT = fullMetadata_OFFS[fullMetadata_OFFS$brotherPairID %in% BP &
                                         fullMetadata_OFFS$trtG1G2 %in% c("NE_exposed", "E_exposed"), ]
  metadataBP_CC_CT = fullMetadata_OFFS[fullMetadata_OFFS$brotherPairID %in% BP &
                                         fullMetadata_OFFS$trtG1G2 %in% c("NE_control", "NE_exposed"), ]
  metadataBP_TC_TT = fullMetadata_OFFS[fullMetadata_OFFS$brotherPairID %in% BP &
                                         fullMetadata_OFFS$trtG1G2 %in% c("E_control", "E_exposed"), ]
  
  ## Make 4 separate uniteCov:
  myuniteCovBP_CC_TC = methylKit::reorganize(methylObj = uniteCovHALF_G2_woSexAndUnknowChrOVERLAP,
                                  treatment = metadataBP_CC_TC$trtG1G2_NUM, sample.ids = metadataBP_CC_TC$ID)
  myuniteCovBP_CT_TT = methylKit::reorganize(methylObj = uniteCovHALF_G2_woSexAndUnknowChrOVERLAP,
                                  treatment = metadataBP_CT_TT$trtG1G2_NUM, sample.ids = metadataBP_CT_TT$ID)
  myuniteCovBP_CC_CT = methylKit::reorganize(methylObj = uniteCovHALF_G2_woSexAndUnknowChrOVERLAP,
                                  treatment = metadataBP_CC_CT$trtG1G2_NUM, sample.ids = metadataBP_CC_CT$ID)
  myuniteCovBP_TC_TT = methylKit::reorganize(methylObj = uniteCovHALF_G2_woSexAndUnknowChrOVERLAP,
                                  treatment = metadataBP_TC_TT$trtG1G2_NUM, sample.ids = metadataBP_TC_TT$ID)
  
  ## remove bases where NO fish in this BP has a coverage
  myuniteCovBP_CC_TC = methylKit::select(myuniteCovBP_CC_TC, which(!is.na(rowSums(percMethylation(myuniteCovBP_CC_TC)))))
  myuniteCovBP_CT_TT = methylKit::select(myuniteCovBP_CT_TT, which(!is.na(rowSums(percMethylation(myuniteCovBP_CT_TT)))))
  myuniteCovBP_CC_CT = methylKit::select(myuniteCovBP_CC_CT, which(!is.na(rowSums(percMethylation(myuniteCovBP_CC_CT)))))
  myuniteCovBP_TC_TT = methylKit::select(myuniteCovBP_TC_TT, which(!is.na(rowSums(percMethylation(myuniteCovBP_TC_TT)))))
  
  ## Calculate differential methylation:
  ## We select the bases that have q-value<0.01 and percent methylation difference larger than 15%, sex as covariate
  DMS_15pc_BP_CC_TC = getDiffMeth(myuniteCov = myuniteCovBP_CC_TC, myMetadata = metadataBP_CC_TC, mccores = 10, mydif = 15)
  DMS_15pc_BP_CT_TT = getDiffMeth(myuniteCov = myuniteCovBP_CT_TT, myMetadata = metadataBP_CT_TT, mccores = 10, mydif = 15)
  DMS_15pc_BP_CC_CT = getDiffMeth(myuniteCov = myuniteCovBP_CC_CT, myMetadata = metadataBP_CC_CT, mccores = 10, mydif = 15)
  DMS_15pc_BP_TC_TT = getDiffMeth(myuniteCov = myuniteCovBP_TC_TT, myMetadata = metadataBP_TC_TT, mccores = 10, mydif = 15)
  
  ## tile for Differentially methylated REGIONS
  tilesBP_CC_TC = methylKit::tileMethylCounts(myuniteCovBP_CC_TC,win.size=100,step.size=100,cov.bases = 10)
  DMR_15pc_BP_CC_TC = getDiffMeth(tilesBP_CC_TC, metadataBP_CC_TC)
  tilesBP_CT_TT = methylKit::tileMethylCounts(myuniteCovBP_CT_TT,win.size=100,step.size=100,cov.bases = 10)
  DMR_15pc_BP_CT_TT = getDiffMeth(tilesBP_CT_TT, metadataBP_CT_TT)
  tilesBP_CC_CT = methylKit::tileMethylCounts(myuniteCovBP_CC_CT,win.size=100,step.size=100,cov.bases = 10)
  DMR_15pc_BP_CC_CT = getDiffMeth(tilesBP_CC_CT, metadataBP_CC_CT)
  tilesBP_TC_TT = methylKit::tileMethylCounts(myuniteCovBP_TC_TT,win.size=100,step.size=100,cov.bases = 10)
  DMR_15pc_BP_TC_TT = getDiffMeth(tilesBP_TC_TT, metadataBP_TC_TT)
  
  return(list(DMSlist = list(DMS_15pc_BP_CC_TC = DMS_15pc_BP_CC_TC, DMS_15pc_BP_CT_TT = DMS_15pc_BP_CT_TT, DMS_15pc_BP_CC_CT = DMS_15pc_BP_CC_CT, DMS_15pc_BP_TC_TT = DMS_15pc_BP_TC_TT),
              DMRlist = list(DMR_15pc_BP_CC_TC = DMR_15pc_BP_CC_TC, DMR_15pc_BP_CT_TT = DMR_15pc_BP_CT_TT, DMR_15pc_BP_CC_CT = DMR_15pc_BP_CC_CT, DMR_15pc_BP_TC_TT = DMR_15pc_BP_TC_TT)))
}

run = FALSE

## Calculate DMS by brother pairs:
if (run == TRUE){

  #### We will apply the following function to all BP:
  vecBP <- unique(fullMetadata_OFFS$brotherPairID)
  
  ## Loop over all BP
  DMlist <- list() # empty plot list
  for (i in 1:length(vecBP)){
    DMlist[[i]] <- getDMperBP(BP = vecBP[[i]])
  } 
  names(DMlist) <- vecBP
  
  saveRDS(DMlist, "../../data/DiffMeth/DMperBP_list_dec2022.RDS")
}

# ## And for positions covered in ALL FISH (not updated)
# getDMperBP2 <- function(BP){
#   ## Unite object for one Brother Pair:
#   metadataBP_CC_TC = fullMetadata[fullMetadata$brotherPairID %in% BP &
#                                     fullMetadata$trtG1G2 %in% c("NE_control", "E_control"), ]
#   metadataBP_CT_TT = fullMetadata[fullMetadata$brotherPairID %in% BP &
#                                     fullMetadata$trtG1G2 %in% c("NE_exposed", "E_exposed"), ]
#   metadataBP_CC_CT = fullMetadata[fullMetadata$brotherPairID %in% BP &
#                                     fullMetadata$trtG1G2 %in% c("NE_control", "NE_exposed"), ]
#   metadataBP_TC_TT = fullMetadata[fullMetadata$brotherPairID %in% BP &
#                                     fullMetadata$trtG1G2 %in% c("E_control", "E_exposed"), ]
#   
#   ## Make 4 separate uniteCov:
#   myuniteCovBP_CC_TC = reorganize(methylObj = uniteCovALL_woSexAndUnknowChr,
#                                   treatment = metadataBP_CC_TC$trtG1G2_NUM, sample.ids = metadataBP_CC_TC$ID)
#   myuniteCovBP_CT_TT = reorganize(methylObj = uniteCovALL_woSexAndUnknowChr,
#                                   treatment = metadataBP_CT_TT$trtG1G2_NUM, sample.ids = metadataBP_CT_TT$ID)
#   myuniteCovBP_CC_CT = reorganize(methylObj = uniteCovALL_woSexAndUnknowChr,
#                                   treatment = metadataBP_CC_CT$trtG1G2_NUM, sample.ids = metadataBP_CC_CT$ID)
#   myuniteCovBP_TC_TT = reorganize(methylObj = uniteCovALL_woSexAndUnknowChr,
#                                   treatment = metadataBP_TC_TT$trtG1G2_NUM, sample.ids = metadataBP_TC_TT$ID)
#   
#   ## remove bases where NO fish in this BP has a coverage
#   myuniteCovBP_CC_TC = methylKit::select(myuniteCovBP_CC_TC, which(!is.na(rowSums(percMethylation(myuniteCovBP_CC_TC)))))
#   myuniteCovBP_CT_TT = methylKit::select(myuniteCovBP_CT_TT, which(!is.na(rowSums(percMethylation(myuniteCovBP_CT_TT)))))
#   myuniteCovBP_CC_CT = methylKit::select(myuniteCovBP_CC_CT, which(!is.na(rowSums(percMethylation(myuniteCovBP_CC_CT)))))
#   myuniteCovBP_TC_TT = methylKit::select(myuniteCovBP_TC_TT, which(!is.na(rowSums(percMethylation(myuniteCovBP_TC_TT)))))
#   
#   ## Calculate differential methylation:
#   ## We select the bases that have q-value<0.01 and percent methylation difference larger than 15%, sex as covariate
#   DMS_15pc_BP_CC_TC = getDiffMeth(myuniteCov = myuniteCovBP_CC_TC, myMetadata = metadataBP_CC_TC, mccores = 10, mydif = 15)
#   DMS_15pc_BP_CT_TT = getDiffMeth(myuniteCov = myuniteCovBP_CT_TT, myMetadata = metadataBP_CT_TT, mccores = 10, mydif = 15)
#   DMS_15pc_BP_CC_CT = getDiffMeth(myuniteCov = myuniteCovBP_CC_CT, myMetadata = metadataBP_CC_CT, mccores = 10, mydif = 15)
#   DMS_15pc_BP_TC_TT = getDiffMeth(myuniteCov = myuniteCovBP_TC_TT, myMetadata = metadataBP_TC_TT, mccores = 10, mydif = 15)
#   
#   ## tile for Differentially methylated REGIONS
#   tilesBP_CC_TC = tileMethylCounts(myuniteCovBP_CC_TC,win.size=100,step.size=100,cov.bases = 10)
#   DMR_15pc_BP_CC_TC = getDiffMeth(tilesBP_CC_TC, metadataBP_CC_TC)
#   tilesBP_CT_TT = tileMethylCounts(myuniteCovBP_CT_TT,win.size=100,step.size=100,cov.bases = 10)
#   DMR_15pc_BP_CT_TT = getDiffMeth(tilesBP_CT_TT, metadataBP_CT_TT)
#   tilesBP_CC_CT = tileMethylCounts(myuniteCovBP_CC_CT,win.size=100,step.size=100,cov.bases = 10)
#   DMR_15pc_BP_CC_CT = getDiffMeth(tilesBP_CC_CT, metadataBP_CC_CT)
#   tilesBP_TC_TT = tileMethylCounts(myuniteCovBP_TC_TT,win.size=100,step.size=100,cov.bases = 10)
#   DMR_15pc_BP_TC_TT = getDiffMeth(tilesBP_TC_TT, metadataBP_TC_TT)
#   
#   return(list(DMSlist = list(DMS_15pc_BP_CC_TC = DMS_15pc_BP_CC_TC, DMS_15pc_BP_CT_TT = DMS_15pc_BP_CT_TT, DMS_15pc_BP_CC_CT = DMS_15pc_BP_CC_CT, DMS_15pc_BP_TC_TT = DMS_15pc_BP_TC_TT),
#               DMRlist = list(DMR_15pc_BP_CC_TC = DMR_15pc_BP_CC_TC, DMR_15pc_BP_CT_TT = DMR_15pc_BP_CT_TT, DMR_15pc_BP_CC_CT = DMR_15pc_BP_CC_CT, DMR_15pc_BP_TC_TT = DMR_15pc_BP_TC_TT)))
# }
# 
# run = FALSE
# if (run == TRUE){
#   #### We will apply the following function to all BP:
#   vecBP <- unique(fullMetadata_OFFS$brotherPairID)
#   
#   ## Loop over all BP
#   DMlist <- list() # empty plot list
#   for (i in 1:length(vecBP)){
#     DMlist[[i]] <- getDMperBP2(BP = vecBP[[i]])
#   } 
#   names(DMlist) <- vecBP
#   
#   saveRDS(DMlist, "../../data/DiffMeth/DMperBP_ALLpos_list.RDS")
# }
