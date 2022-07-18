## Differential methylation analyses
## A. Balard
## February 2022

machine="apocrita" # define the machine we work on
loadALL = FALSE # only load CpG shared by half fish per trt group
loadannot = FALSE # we don't need annotations at this stage
sourceDMS = FALSE # we calculate DMS here
sourceSubUnite = TRUE # get the subunite objects to calculate differential methylations
source("R02.3_DATALOAD.R")

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
run = FALSE

if (run == TRUE){
  DMS_15pc_G1_C_T = getDiffMeth(uniteCov6_G1_woSexAndUnknowChrOVERLAP, fullMetadata_PAR) # old name: DMS15pc_G1_half
  DMS_15pc_CC_TC = getDiffMeth(uniteCov14_G1bothTrt_G2control_woSexAndUnknowChr, fullMetadata_OFFS[fullMetadata_OFFS$trtG1G2_NUM %in% c(2,5),]) # old DMS15pc_G1_controlG2_half
  DMS_15pc_CT_TT = getDiffMeth(uniteCov14_G1bothTrt_G2infected_woSexAndUnknowChr, fullMetadata_OFFS[fullMetadata_OFFS$trtG1G2_NUM %in% c(3,6),]) # old DMS15pc_G1_infectedG2_half
  DMS_15pc_CC_CT = getDiffMeth(uniteCov14_G2_woSexAndUnknowChr_G1CONTROL, fullMetadata_OFFS[fullMetadata_OFFS$trtG1G2_NUM %in% c(5,6),]) # old name DMS15pc_G2_controlG1_half
  DMS_15pc_TC_TT = getDiffMeth(uniteCov14_G2_woSexAndUnknowChr_G1INFECTED, fullMetadata_OFFS[fullMetadata_OFFS$trtG1G2_NUM %in% c(2,3),]) # old DMS15pc_G2_infectedG1_half
  
  DMR_15pc_G1_C_T = getDiffMeth(tiles_G1_half, fullMetadata_PAR) # old name: DMR15pc_G1_half
  DMR_15pc_CC_TC = getDiffMeth(tiles_G1bothTrt_G2control_half, fullMetadata_OFFS[fullMetadata_OFFS$trtG1G2_NUM %in% c(2,5),])
  DMR_15pc_CT_TT = getDiffMeth(tiles_G1bothTrt_G2infected_half, fullMetadata_OFFS[fullMetadata_OFFS$trtG1G2_NUM %in% c(3,6),])
  DMR_15pc_CC_CT = getDiffMeth(tiles_G2_G1CONTROL_half, fullMetadata_OFFS[fullMetadata_OFFS$trtG1G2_NUM %in% c(5,6),])
  DMR_15pc_TC_TT = getDiffMeth(tiles_G2_G1INFECTED_half, fullMetadata_OFFS[fullMetadata_OFFS$trtG1G2_NUM %in% c(2,3),])
  
  # Save list of DMS/DMR
  DMSlist = list(DMS_15pc_G1_C_T = DMS_15pc_G1_C_T, 
                 DMS_15pc_CC_TC = DMS_15pc_CC_TC, DMS_15pc_CT_TT = DMS_15pc_CT_TT, 
                 DMS_15pc_CC_CT = DMS_15pc_CC_CT, DMS_15pc_TC_TT = DMS_15pc_TC_TT)
  saveRDS(DMSlist, "../../data/DiffMeth/DMSlist.RDS")
  
  DMRlist = list(DMR_15pc_G1_C_T = DMR_15pc_G1_C_T, 
                 DMR_15pc_CC_TC = DMR_15pc_CC_TC, DMR_15pc_CT_TT = DMR_15pc_CT_TT, 
                 DMR_15pc_CC_CT = DMR_15pc_CC_CT, DMR_15pc_TC_TT = DMR_15pc_TC_TT)
  saveRDS(DMRlist, "../../data/DiffMeth/DMRlist.RDS")
}

######################## 
## And by brother pairs:

run = FALSE
if (run == TRUE){
  #### We will apply the following function to all BP:
  vecBP <- unique(fullMetadata_OFFS$brotherPairID)
  
  ## Loop over all BP
  DMlist <- list() # empty plot list
  for (i in 1:length(vecBP)){
    DMlist[[i]] <- getDMperBP(BP = vecBP[[i]])
  } 
  names(DMlist) <- vecBP
  
  saveRDS(DMlist, "../../data/DiffMeth/DMperBP_list.RDS")
}

## For positions present in all fish (for upset plot)
run = TRUE
if (run == TRUE){
  #### We will apply the following function to all BP:
  vecBP <- unique(fullMetadata_OFFS$brotherPairID)
  
  ## Loop over all BP
  DMlist <- list() # empty plot list
  for (i in 1:length(vecBP)){
    DMlist[[i]] <- getDMperBP2(BP = vecBP[[i]])
  } 
  names(DMlist) <- vecBP
  
  saveRDS(DMlist, "../../data/DiffMeth/DMperBP_ALLpos_list.RDS")
}
