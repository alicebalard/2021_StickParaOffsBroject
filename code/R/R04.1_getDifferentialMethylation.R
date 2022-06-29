## Differential methylation analyses
## A. Balard
## February 2022

machine="apocrita" # define the machine we work on
loadALL = FALSE # only load CpG shared by half fish per trt group
loadannot = FALSE # we don't need annotations at this stage
sourceDMS = FALSE # we calculate DMS here
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

######################################
## For DMR: tile the methylation data:
## Kostas MBE 2020: "To identify DMRs, we used the tileMethylCounts() function in MethylKit
# v.1.5.0 with a sliding window size of 100 bases and step size of 100 bases."
# Summarize methylation information over tiling windows with a sliding window size of 100 bases and step size of 100 bases

# 0. PARENTS trt-ctrl (CpG covered in half trt group)
tiles_G1_half = tileMethylCounts(uniteCov6_G1_woSexAndUnknowChrOVERLAP, win.size=100,step.size=100,cov.bases = 10)
nrow(tiles_G1_half) # methylBase object with 20348 rows

# 1. CC-TC = CONTROL fish (parent CvsT) (CpG covered in half trt group)
uniteCov14_G1bothTrt_G2control_woSexAndUnknowChr <- reorganize(methylObj = uniteCov14_G2_woSexAndUnknowChrOVERLAP,
                                                               treatment = fullMetadata_OFFS$trtG1G2_NUM[fullMetadata_OFFS$trtG1G2_NUM %in% c(2,5)],
                                                               sample.ids = fullMetadata_OFFS$ID[fullMetadata_OFFS$trtG1G2_NUM %in% c(2,5)])
tiles_G1bothTrt_G2control_half = tileMethylCounts(uniteCov14_G1bothTrt_G2control_woSexAndUnknowChr,win.size=100,step.size=100,cov.bases = 10)
nrow(tiles_G1bothTrt_G2control_half) # methylBase object with 20348 rows

# 2. CT-TT = TREATMENT fish (parent CvsT) (CpG covered in half trt group)
uniteCov14_G1bothTrt_G2infected_woSexAndUnknowChr <- reorganize(methylObj = uniteCov14_G2_woSexAndUnknowChrOVERLAP,
                                                                treatment = fullMetadata_OFFS$trtG1G2_NUM[fullMetadata_OFFS$trtG1G2_NUM %in% c(3,6)],
                                                                sample.ids = fullMetadata_OFFS$ID[fullMetadata_OFFS$trtG1G2_NUM %in% c(3,6)])
tiles_G1bothTrt_G2infected_half = tileMethylCounts(uniteCov14_G1bothTrt_G2infected_woSexAndUnknowChr,win.size=100,step.size=100,cov.bases = 10)
nrow(tiles_G1bothTrt_G2infected_half) # methylBase object with 20348 rows

# 3. CC-CT = fish from CONTROL parents (G2 CvsT) (CpG covered in half trt group)
uniteCov14_G2_woSexAndUnknowChr_G1CONTROL <- reorganize(methylObj = uniteCov14_G2_woSexAndUnknowChrOVERLAP,
                                                        treatment = fullMetadata_OFFS$trtG1G2_NUM[fullMetadata_OFFS$trtG1G2_NUM %in% c(5,6)],
                                                        sample.ids = fullMetadata_OFFS$ID[fullMetadata_OFFS$trtG1G2_NUM %in% c(5,6)])
tiles_G2_G1CONTROL_half = tileMethylCounts(uniteCov14_G2_woSexAndUnknowChr_G1CONTROL,win.size=100,step.size=100,cov.bases = 10)
nrow(tiles_G2_G1CONTROL_half) # methylBase object with 20348 rows

# 4. TC-TT = fish from TREATMENT parents (G2 CvsT) (CpG covered in half trt group)
uniteCov14_G2_woSexAndUnknowChr_G1INFECTED <- reorganize(methylObj = uniteCov14_G2_woSexAndUnknowChrOVERLAP,
                                                         treatment = fullMetadata_OFFS$trtG1G2_NUM[fullMetadata_OFFS$trtG1G2_NUM %in% c(2,3)],
                                                         sample.ids = fullMetadata_OFFS$ID[fullMetadata_OFFS$trtG1G2_NUM %in% c(2,3)])
tiles_G2_G1INFECTED_half = tileMethylCounts(uniteCov14_G2_woSexAndUnknowChr_G1INFECTED,win.size=100,step.size=100,cov.bases = 10)
nrow(tiles_G2_G1INFECTED_half) # methylBase object with 20348 rows

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
    DMSlist = list(DMS_15pc_G1_C_T, DMS_15pc_CC_TC, DMS_15pc_CT_TT, DMS_15pc_CC_CT, DMS_15pc_TC_TT)
                                        # saveRDS(DMSlist, "../../data/DiffMeth/DMSlist.RDS")

    DMRlist = list(DMR_15pc_G1_C_T, DMR_15pc_CC_TC, DMR_15pc_CT_TT, DMR_15pc_CC_CT, DMR_15pc_TC_TT)
                                        # saveRDS(DMRlist, "../../data/DiffMeth/DMRlist.RDS")
}

######################## 
## And by brother pairs:

#### We will apply the following function to all BP:
vecBP <- unique(fullMetadata_OFFS$brotherPairID)

## Loop over all BP
DMlist <- list() # empty plot list
for (i in 1:length(vecBP)){
  DMlist[[i]] <- getDMperBP(BP = vecBP[[i]])
} 
names(DMlist) <- vecBP

saveRDS(DMlist, "../../data/DiffMeth/DMperBP_list.RDS")






#############################
### RM later if found useless
# ## CpG covered in all fish, parental DMS:
# uniteCovALL_G1_woSexAndUnknowChr <- reorganize(methylObj = uniteCovALL_woSexAndUnknowChr,
#                                                treatment = fullMetadata$trtG1G2_NUM[fullMetadata$trtG1G2_NUM %in% c(1,4)],
#                                                sample.ids = fullMetadata$ID[fullMetadata$trtG1G2_NUM %in% c(1,4)])
# nrow(uniteCovALL_G1_woSexAndUnknowChr) # methylBase object with 55530 rows (CpG covered)
# tiles_G1_ALL = tileMethylCounts(uniteCovALL_G1_woSexAndUnknowChr,win.size=100,step.size=100,cov.bases = 10)
# nrow(tiles_G1_ALL) # methylBase object with 413 rows

##############
# ## CpG covered in all fish, G2 from G1 control DMS:
# uniteCovALL_G2_woSexAndUnknowChr_G1CONTROL <- reorganize(methylObj = uniteCovALL_woSexAndUnknowChr,
#                                                          treatment = fullMetadata$trtG1G2_NUM[
#                                                            fullMetadata$trtG1G2_NUM %in% c(5,6)],
#                                                          sample.ids = fullMetadata$ID[
#                                                            fullMetadata$trtG1G2_NUM %in% c(5,6)])
# nrow(uniteCovALL_G2_woSexAndUnknowChr_G1CONTROL) # methylBase object with 55530 rows (CpG covered)
# 
# tiles_G2_G1CONTROL_ALL = tileMethylCounts(uniteCovALL_G2_woSexAndUnknowChr_G1CONTROL,win.size=100,step.size=100,cov.bases = 10)
# nrow(tiles_G2_G1CONTROL_ALL) # methylBase object with 413 rows

##############
# ## CpG covered in all fish, G2 from G1 infected DMS:
# uniteCovALL_G2_woSexAndUnknowChr_G1INFECTED <- reorganize(methylObj = uniteCovALL_woSexAndUnknowChr,
#                                                           treatment = fullMetadata$trtG1G2_NUM[
#                                                             fullMetadata$trtG1G2_NUM %in% c(2,3)],
#                                                           sample.ids = fullMetadata$ID[
#                                                             fullMetadata$trtG1G2_NUM %in% c(2,3)])
# nrow(uniteCovALL_G2_woSexAndUnknowChr_G1INFECTED) # methylBase object with 55530 rows (CpG covered)
# tiles_G2_G1INFECTED_ALL = tileMethylCounts(uniteCovALL_G2_woSexAndUnknowChr_G1INFECTED,win.size=100,step.size=100,cov.bases = 10)
# nrow(tiles_G2_G1INFECTED_ALL) # methylBase object with 413 rows

## stop here:
## stop("We stop here for now") # to run getDiffMeth on Apocrita cause it's LONG

## To do on bash: rename with DATES 
## for f in DMS15pc*; do mv "$f" "$(echo "$f" | sed s/.RDS/_25feb22.RDS/)"; done
## for f in DMR15pc*; do mv "$f" "$(echo "$f" | sed s/.RDS/_25feb22.RDS/)"; done

