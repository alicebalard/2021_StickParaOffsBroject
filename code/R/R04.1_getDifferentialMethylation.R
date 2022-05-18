## Differential methylation analyses
## A. Balard
## February 2022

machine="apocrita" # define the machine we work on
loadALL = FALSE # only load CpG shared by half fish per trt group
source("R02.3_DATALOAD.R")

#######################################################
### PART 2: Differential Methylation Sites/Regions ####
#######################################################

### 3 comparisons to do:
## PARENTS trt-ctrl
## G1control G2 trt-ctrl
## G1infected G2 trt-ctrl ## HYP: this will be different

## Extra 2 for paternal effect:
## G1control-G2control vs G1infected-G2control
## G1control-G2infected vs G1infected-G2infected

################################################################
## Correspondance trtG1G2 and numerical values used by MethylKit
#table(fullMetadata$trtG1G2_NUM, fullMetadata$trtG1G2)
# Control Exposed NE_control NE_exposed E_control E_exposed
# 1      12       0          0          0         0         0
# 2       0       0          0          0        28         0
# 3       0       0          0          0         0        28
# 4       0      12          0          0         0         0
# 5       0       0         28          0         0         0
# 6       0       0          0         27         0         0

######################################
## For DMR: tile the methylation data:
## Kostas MBE 2020: "To identify DMRs, we used the tileMethylCounts() function in MethylKit
# v.1.5.0 with a sliding window size of 100 bases and step size of 100 bases."

## Summarize methylation information over tiling windows with a sliding window size of 100 bases and step size of 100 bases

################
##### G1
##### CpG covered in half trt group, parental DMS:
###nrow(uniteCov6_G1_woSexAndUnknowChrOVERLAP) # methylBase object with 1001880 rows
###tiles_G1_half = tileMethylCounts(uniteCov6_G1_woSexAndUnknowChrOVERLAP, win.size=100,step.size=100,cov.bases = 10)
###nrow(tiles_G1_half) # methylBase object with 20348 rows
###
##### CpG covered in all fish, parental DMS:
###uniteCovALL_G1_woSexAndUnknowChr <- reorganize(methylObj = uniteCovALL_woSexAndUnknowChr,
###                   treatment = fullMetadata$trtG1G2_NUM[
###                                                     fullMetadata$trtG1G2_NUM %in% c(1,4)],
###                   sample.ids = fullMetadata$ID[
###                                                 fullMetadata$trtG1G2_NUM %in% c(1,4)])
###nrow(uniteCovALL_G1_woSexAndUnknowChr) # methylBase object with 55530 rows (CpG covered)
###tiles_G1_ALL = tileMethylCounts(uniteCovALL_G1_woSexAndUnknowChr,win.size=100,step.size=100,cov.bases = 10)
###nrow(tiles_G1_ALL) # methylBase object with 413 rows
###
#################
##### G2 from control fathers
##### CpG covered in half trt group, G2 from G1 control DMS:
###uniteCov14_G2_woSexAndUnknowChr_G1CONTROL <- reorganize(methylObj = uniteCov14_G2_woSexAndUnknowChrOVERLAP,
###                                                        treatment = fullMetadata_OFFS$trtG1G2_NUM[
###                                                          fullMetadata_OFFS$trtG1G2_NUM %in% c(5,6)],
###                                                        sample.ids = fullMetadata_OFFS$ID[
###                                                          fullMetadata_OFFS$trtG1G2_NUM %in% c(5,6)])
###nrow(uniteCov14_G2_woSexAndUnknowChr_G1CONTROL) # methylBase object with 1001880 rows
###tiles_G2_G1CONTROL_half = tileMethylCounts(uniteCov14_G2_woSexAndUnknowChr_G1CONTROL,win.size=100,step.size=100,cov.bases = 10)
###nrow(tiles_G2_G1CONTROL_half) # methylBase object with 20348 rows
###
##### CpG covered in all fish, G2 from G1 control DMS:
###uniteCovALL_G2_woSexAndUnknowChr_G1CONTROL <- reorganize(methylObj = uniteCovALL_woSexAndUnknowChr,
###                   treatment = fullMetadata$trtG1G2_NUM[
###                                                     fullMetadata$trtG1G2_NUM %in% c(5,6)],
###                   sample.ids = fullMetadata$ID[
###                                                 fullMetadata$trtG1G2_NUM %in% c(5,6)])
###nrow(uniteCovALL_G2_woSexAndUnknowChr_G1CONTROL) # methylBase object with 55530 rows (CpG covered)
###
###tiles_G2_G1CONTROL_ALL = tileMethylCounts(uniteCovALL_G2_woSexAndUnknowChr_G1CONTROL,win.size=100,step.size=100,cov.bases = 10)
###nrow(tiles_G2_G1CONTROL_ALL) # methylBase object with 413 rows
###
#################
##### G2 from infected fathers
##### CpG covered in half trt group, G2 from G1 infected DMS:
###uniteCov14_G2_woSexAndUnknowChr_G1INFECTED <- reorganize(methylObj = uniteCov14_G2_woSexAndUnknowChrOVERLAP,
###                                                        treatment = fullMetadata_OFFS$trtG1G2_NUM[
###                                                          fullMetadata_OFFS$trtG1G2_NUM %in% c(2,3)],
###                                                        sample.ids = fullMetadata_OFFS$ID[
###                                                          fullMetadata_OFFS$trtG1G2_NUM %in% c(2,3)])
###nrow(uniteCov14_G2_woSexAndUnknowChr_G1INFECTED) # methylBase object with 1001880 rows
###
###tiles_G2_G1INFECTED_half = tileMethylCounts(uniteCov14_G2_woSexAndUnknowChr_G1INFECTED,win.size=100,step.size=100,cov.bases = 10)
###nrow(tiles_G2_G1INFECTED_half) # methylBase object with 20348 rows
###
##### CpG covered in all fish, G2 from G1 infected DMS:
###uniteCovALL_G2_woSexAndUnknowChr_G1INFECTED <- reorganize(methylObj = uniteCovALL_woSexAndUnknowChr,
###                   treatment = fullMetadata$trtG1G2_NUM[
###                                                     fullMetadata$trtG1G2_NUM %in% c(2,3)],
###                   sample.ids = fullMetadata$ID[
###                                                 fullMetadata$trtG1G2_NUM %in% c(2,3)])
###nrow(uniteCovALL_G2_woSexAndUnknowChr_G1INFECTED) # methylBase object with 55530 rows (CpG covered)
###tiles_G2_G1INFECTED_ALL = tileMethylCounts(uniteCovALL_G2_woSexAndUnknowChr_G1INFECTED,win.size=100,step.size=100,cov.bases = 10)
###nrow(tiles_G2_G1INFECTED_ALL) # methylBase object with 413 rows
###

##############
## G2 control from both control & infected fathers
## CpG covered in half trt group
uniteCov14_G1bothTrt_G2control_woSexAndUnknowChr <- reorganize(methylObj = uniteCov14_G2_woSexAndUnknowChrOVERLAP,
                                                        treatment = fullMetadata_OFFS$trtG1G2_NUM[
                                                          fullMetadata_OFFS$trtG1G2_NUM %in% c(2,5)],
                                                        sample.ids = fullMetadata_OFFS$ID[
                                                          fullMetadata_OFFS$trtG1G2_NUM %in% c(2,5)])
nrow(uniteCov14_G1bothTrt_G2control_woSexAndUnknowChr) # methylBase object with 1001880 rows

tiles_G1bothTrt_G2control_half = tileMethylCounts(uniteCov14_G1bothTrt_G2control_woSexAndUnknowChr,win.size=100,step.size=100,cov.bases = 10)
nrow(tiles_G1bothTrt_G2control_half) # methylBase object with 20348 rows

## G2 infected from both control & infected fathers
## CpG covered in half trt group
uniteCov14_G1bothTrt_G2infected_woSexAndUnknowChr <- reorganize(methylObj = uniteCov14_G2_woSexAndUnknowChrOVERLAP,
                                                        treatment = fullMetadata_OFFS$trtG1G2_NUM[
                                                          fullMetadata_OFFS$trtG1G2_NUM %in% c(3,6)],
                                                        sample.ids = fullMetadata_OFFS$ID[
                                                          fullMetadata_OFFS$trtG1G2_NUM %in% c(3,6)])
nrow(uniteCov14_G1bothTrt_G2infected_woSexAndUnknowChr) # methylBase object with 1001880 rows

tiles_G1bothTrt_G2infected_half = tileMethylCounts(uniteCov14_G1bothTrt_G2infected_woSexAndUnknowChr,win.size=100,step.size=100,cov.bases = 10)
nrow(tiles_G1bothTrt_G2infected_half) # methylBase object with 20348 rows

##############################################################################
## Calculate DMS/DMR accounting for covariates: brotherPairID and sex (new 04/04/22!)

###### Comparison 1: BASELINE -> Parents (control vs infected) 
###### CpG covered in HALF fish per group
####DMS15pc_G1_half <- getDiffMeth(uniteCov6_G1_woSexAndUnknowChrOVERLAP, fullMetadata_PAR)
####saveRDS(DMS15pc_G1_half, file = "../../data/DiffMeth/DMS15pc_G1_half.RDS")
####
####DMR15pc_G1_half <- getDiffMeth(tiles_G1_half, fullMetadata_PAR)
####saveRDS(DMR15pc_G1_half, file = "../../data/DiffMeth/DMR15pc_G1_half.RDS")
####
###### CpG covered in ALL fish
####DMS15pc_G1_ALL <- getDiffMeth(uniteCovALL_G1_woSexAndUnknowChr, fullMetadata_PAR)
####saveRDS(DMS15pc_G1_ALL, file = "../../data/DiffMeth/DMS15pc_G1_ALL.RDS")
####
####DMR15pc_G1_ALL <- getDiffMeth(tiles_G1_ALL, fullMetadata_PAR)
####DMR15pc_G1_ALL # EMPTY
####
###### Comparison 2: Offspring
###### 2.1. Should be like baseline -> G2 from G1 control (control vs infected)
####
###### CpG covered in HALF fish per group
####DMS15pc_G2_controlG1_half <- getDiffMeth(myuniteCov = uniteCov14_G2_woSexAndUnknowChr_G1CONTROL,
####                                     myMetadata = fullMetadata_OFFS[fullMetadata_OFFS$trtG1G2_NUM %in% c(5,6),])
####saveRDS(DMS15pc_G2_controlG1_half, file = "../../data/DiffMeth/DMS15pc_G2_controlG1_half.RDS")
####
####DMR15pc_G2_controlG1_half <- getDiffMeth(tiles_G2_G1CONTROL_half, fullMetadata_OFFS[fullMetadata_OFFS$trtG1G2_NUM %in% c(5,6),])
####saveRDS(DMR15pc_G2_controlG1_half, file = "../../data/DiffMeth/DMR15pc_G2_controlG1_half.RDS")
####
###### CpG covered in ALL fish
####DMS15pc_G2_controlG1_ALL <- getDiffMeth(myuniteCov = uniteCovALL_G2_woSexAndUnknowChr_G1CONTROL,
####                                        myMetadata = fullMetadata_OFFS[fullMetadata_OFFS$trtG1G2_NUM %in% c(5,6),])
####
####DMS15pc_G2_controlG1_ALL
####saveRDS(DMS15pc_G2_controlG1_ALL, file = "../../data/DiffMeth/DMS15pc_G2_controlG1_ALL.RDS")
####
####DMR15pc_G2_controlG1_ALL <- getDiffMeth(tiles_G2_G1CONTROL_ALL, fullMetadata_OFFS[fullMetadata_OFFS$trtG1G2_NUM %in% c(5,6),])
####saveRDS(DMR15pc_G2_controlG1_ALL, file = "../../data/DiffMeth/DMR15pc_G2_controlG1_ALL.RDS")
####
###### 2.2. Should be DIFFERENT if there is a paternal effect -> G2 from G1 infected (control vs infected) 
####
###### CpG covered in HALF fish per group
####DMS15pc_G2_infectedG1_half <- getDiffMeth(myuniteCov = uniteCov14_G2_woSexAndUnknowChr_G1INFECTED, 
####                                      myMetadata = fullMetadata_OFFS[fullMetadata_OFFS$trtG1G2_NUM %in% c(2,3),])
####  saveRDS(DMS15pc_G2_infectedG1_half, file = "../../data/DiffMeth/DMS15pc_G2_infectedG1_half.RDS")
####
####DMR15pc_G2_infectedG1_half <- getDiffMeth(tiles_G2_G1INFECTED_half, fullMetadata_OFFS[fullMetadata_OFFS$trtG1G2_NUM %in% c(2,3),])
####saveRDS(DMR15pc_G2_infectedG1_half, file = "../../data/DiffMeth/DMR15pc_G2_infectedG1_half.RDS")
####
###### CpG covered in ALL fish
####DMS15pc_G2_infectedG1_ALL <- getDiffMeth(myuniteCov = uniteCovALL_G2_woSexAndUnknowChr_G1INFECTED,
####                                        myMetadata = fullMetadata_OFFS[fullMetadata_OFFS$trtG1G2_NUM %in% c(2,3),])
####saveRDS(DMS15pc_G2_infectedG1_ALL, file = "../../data/DiffMeth/DMS15pc_G2_infectedG1_ALL.RDS")
####
####DMR15pc_G2_controlG1_ALL <- getDiffMeth(tiles_G2_G1INFECTED_ALL, fullMetadata_OFFS[fullMetadata_OFFS$trtG1G2_NUM %in% c(2,3),])
####saveRDS(DMR15pc_G2_controlG1_ALL, file = "../../data/DiffMeth/DMR15pc_G2_infectedG1_ALL.RDS")
####

### 3.1. Should be DIFFERENT if there is a paternal effect -> G2 control from both G1; CpG covered in HALF fish per group
## DMS:
DMS15pc_G1_controlG2_half <- getDiffMeth(myuniteCov = uniteCov14_G1bothTrt_G2control_woSexAndUnknowChr,
                                      myMetadata = fullMetadata_OFFS[fullMetadata_OFFS$trtG1G2_NUM %in% c(2,5),])
saveRDS(DMS15pc_G1_controlG2_half, file = "../../data/DiffMeth/DMS15pc_G1_controlG2_half.RDS")

## DMR:
DMR15pc_G1_controlG2_half <- getDiffMeth(tiles_G1bothTrt_G2control_half, fullMetadata_OFFS[fullMetadata_OFFS$trtG1G2_NUM %in% c(2,5),])
saveRDS(DMR15pc_G1_controlG2_half, file = "../../data/DiffMeth/DMR15pc_G1_controlG2_half.RDS")

### 3.2. Should be DIFFERENT if there is a paternal effect -> G2 infected from both G1; CpG covered in HALF fish per group
##DMS:
DMS15pc_G1_infectedG2_half <- getDiffMeth(myuniteCov = uniteCov14_G1bothTrt_G2infected_woSexAndUnknowChr,
                                      myMetadata = fullMetadata_OFFS[fullMetadata_OFFS$trtG1G2_NUM %in% c(3,6),])
saveRDS(DMS15pc_G1_infectedG2_half, file = "../../data/DiffMeth/DMS15pc_G1_infectedG2_half.RDS")

## DMR:
DMR15pc_G1_infectedG2_half <- getDiffMeth(tiles_G1bothTrt_G2infected_half, fullMetadata_OFFS[fullMetadata_OFFS$trtG1G2_NUM %in% c(3,6),])
saveRDS(DMR15pc_G1_infectedG2_half, file = "../../data/DiffMeth/DMR15pc_G1_infectedG2_half.RDS")

## stop here:
## stop("We stop here for now") # to run getDiffMeth on Apocrita cause it's LONG

## To do on bash: rename with DATES 
## for f in DMS15pc*; do mv "$f" "$(echo "$f" | sed s/.RDS/_25feb22.RDS/)"; done
## for f in DMR15pc*; do mv "$f" "$(echo "$f" | sed s/.RDS/_25feb22.RDS/)"; done
