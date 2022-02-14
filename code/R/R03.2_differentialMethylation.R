## Differential methylation analyses
## A. Balard
## February 2022

#################### Data load & preparation ####################
source("librariesLoading.R")
## load custom functions
source("customRfunctions.R")
## Load samples metadata
source("R01.3_loadMetadata.R")
## define in which machine we're working (apocrita or mythinkpad)
# machine="apocrita"
machine="mythinkpad"
## Load methylation data
source("R01.4_loadMethyldata.R")

###############################################
### PART 2: Differential methylation sites ####
###############################################

### 3 comparisons to do:
## PARENTS trt-ctrl
## G2-control G1 trt-ctrl
## G2-infected G1 trt-ctrl ## HYP: this will be different

## Calculate DMS accounting for covariates: family
getDMS <- function(myuniteCov, myMetadata){
  cov = data.frame(Family = myMetadata$Family)
  myDiffMeth=calculateDiffMeth(myuniteCov, covariates = cov, mc.cores = 10)
  # We select the bases that have q-value<0.01 and percent methylation difference larger than 15%.
  # NB: arg type="hyper" or type="hypo" gives hyper-methylated or hypo-methylated regions/bases.
  myDMS_15pc = getMethylDiff(myDiffMeth, difference=15, qvalue=0.01)
  return(myDMS_15pc)
}

##########################################################
## Comparison 1: BASELINE -> Parents (control vs infected) 
## Extract DMS from parents (at least in half the fish), annotate them, compare with Kostas results

## rerun or upload
### RERUN
# DMS15pc_PAR_half <- getDMS(uniteCov6_G1_woSexAndUnknowChr, fullMetadata_PAR_half)
# saveRDS(DMS15pc_PAR_half, file = "../../data/DMS15pc_PAR_half.RDS")
### UPLOAD
DMS15pc_PAR_half <- readRDS("../../data/DMS15pc_PAR_half.RDS")
nrow(DMS15pc_PAR_half)
# 6603 positions

## NB Kostas' results: "We found a total of 1,973 CpG sites out of
# 1,172,887 CpGs (0.17%) across the genome that showed at
# least 15% differential fractional methylation (differentially
# methylated site [DMS]; q < 0.01) between infected and uninfected fish"

nrow(uniteCov6_G1_woSexAndUnknowChr)#1 188 179 CpG sites
nrow(DMS15pc_PAR_half) # 6603 DMS
nrow(DMS15pc_PAR_half) / nrow(uniteCov6_G1_woSexAndUnknowChr) *100 # 0.55%

##########################################################
## Comparison 2: Should be like baseline -> G2 from G1 control (control vs infected) 

table(fullMetadata_OFFS$trtG1G2, fullMetadata_OFFS$trtG1G2_NUM)
#             2  3  5  6
# NE_control  0  0 28  0
# NE_exposed  0  0  0 27
# E_control  28  0  0  0
# E_exposed   0 28  0  0

# DMS15pc_G2_controlG1_half <- getDMS(myuniteCov = reorganize(methylObj = uniteCov14_G2_woSexAndUnknowChr,
#                                                             treatment = fullMetadata_OFFS$trtG1G2_NUM[
#                                                                                               fullMetadata_OFFS$trtG1G2_NUM %in% c(5,6)], 
#                                                             sample.ids = fullMetadata_OFFS$ID[
#                                                                                                fullMetadata_OFFS$trtG1G2_NUM %in% c(5,6)]), 
#                                     myMetadata = fullMetadata_OFFS[fullMetadata_OFFS$trtG1G2_NUM %in% c(5,6),])
# 
# saveRDS(DMS15pc_G2_controlG1_half, file = "../../data/DMS15pc_G2_controlG1_half.RDS")

# DMS15pc_G2_infectedG1_half <- getDMS(myuniteCov = reorganize(methylObj = uniteCov14_G2_woSexAndUnknowChr,
#                                                             treatment = fullMetadata_OFFS$trtG1G2_NUM[
#                                                                                               fullMetadata_OFFS$trtG1G2_NUM %in% c(2,3)], 
#                                                             sample.ids = fullMetadata_OFFS$ID[
#                                                                                                fullMetadata_OFFS$trtG1G2_NUM %in% c(2,3)]), 
#                                     myMetadata = fullMetadata_OFFS[fullMetadata_OFFS$trtG1G2_NUM %in% c(2,3),])
# 
# saveRDS(DMS15pc_G2_infectedG1_half, file = "../../data/DMS15pc_G2_infectedG1_half.RDS")

## stop here:
# stop("We stop here for now") # to run getDMS on Apocrita cause it's LONG

### UPLOAD
## Control G1 - G2(trt vs control)
DMS15pc_G2_controlG1_half <- readRDS("../../data/DMS15pc_G2_controlG1_half.RDS")
nrow(DMS15pc_G2_controlG1_half)
# 1642 positions

## Infected G1 - G2(trt vs control)
DMS15pc_G2_infectedG1_half <- readRDS("../../data/DMS15pc_G2_infectedG1_half.RDS")
nrow(DMS15pc_G2_infectedG1_half)
# 943 positions

#### WHICH positions? 
library(VennDiagram)#to put up

## Check which sequenced CpG are overlapping between half offsprings and half parents datasets:

## YOU'RE HERE 14th feb evening!!


## Check which DMS are averlapping between offsprings and parents datasets

posG1DMS <- paste(DMS15pc_PAR_half$chr, DMS15pc_PAR_half$start, DMS15pc_PAR_half$end)
posG2DMS_controlG1 <- paste(DMS15pc_G2_controlG1_half$chr, DMS15pc_G2_controlG1_half$start, DMS15pc_G2_controlG1_half$end)
posG2DMS_infectedG1 <- paste(DMS15pc_G2_infectedG1_half$chr, DMS15pc_G2_infectedG1_half$start, DMS15pc_G2_infectedG1_half$end)

VennDiagram::get.venn.partitions(
  list(posG1DMS=posG1DMS, posG2DMS_controlG1=posG2DMS_controlG1, posG2DMS_infectedG1=posG2DMS_infectedG1))
grid.newpage()
grid::grid.draw(VennDiagram::venn.diagram(list(posG1DMS=posG1DMS, posG2DMS_controlG1=posG2DMS_controlG1, posG2DMS_infectedG1=posG2DMS_infectedG1), NULL))



#############
## annotation

## load genome annotation
gene.obj=readTranscriptFeatures("../../gitignore/bigdata/Gy_allnoM_rd3.maker_apocrita.noseq_corrected.bed12", 
                                remove.unusual = FALSE)
 
## PARENTS trt-ctrl
## G2-control G1 trt-ctrl
## G2-infected G1 trt-ctrl

# annotate differentially methylated CpGs with promoter/exon/intron using annotation data
diffAnn_PAR=annotateWithGeneParts(as(DMS15pc_PAR_half,"GRanges"),gene.obj)
diffAnn_PAR
plotTargetAnnotation(diffAnn_PAR,precedence=TRUE,
                     main="differential methylation annotation")
## percentage of target features overlapping with annotation:
# (with promoter > exon > intron precedence):
# promoter  exon       intron   intergenic 
# 9.13      12.75      30.80      47.31 

## Kostas MBE: The DMSs and regions were predominately
# found in intergenic regions (47.74% and 48.94%, respectively),
# with introns (26.19% and 23.09), exons (15.07% and 13.98%), and promoters (11% and 13.98%) showing lower proportions
 
###########################
## Manhattan plot of DMS ##
###########################

## Parents trt-ctrl

# load annotation
annot_PAR <- as.data.frame(diffAnn_PAR@members)
## Load file containing length of each gynogen chromosomes 
## grep "contig" gitignore/bigdata/Gy_allnoM_rd3.maker_apocrita.noseq_corrected.gff | awk '{print $1, $5}' > data/Gy_allnoM_rd3.maker_apocrita.noseq_corrected_chromoAndLength.txt
GYgynogff = read.table("../../data/Gy_allnoM_rd3.maker_apocrita.noseq_corrected_chromoAndLength.txt")
names(GYgynogff) = c("chrom","length")

makeManhattanPlots(DMSfile = DMS15pc_PAR_half, annotFile = annot_PAR, GYgynogff = GYgynogff, 
                   mycols = c("red", "grey", "black", "green"))


# GO terms --> TBC
# Networks --> TBC


# 
# ###################################################################################################
# ## Comparison 2: offsprings infected from unifected (trtgroup 4) and infected (trt group 6) fathers
# ## NB CpG present in ALL FISH (not at least 2) for a start
# 
# ### Select metadata
# myMetaData = fullMetadata_OFFS[fullMetadata_OFFS$trtG1G2_NUM %in% c(4,6),]
# 
# ### Select methylKit object
# myMethylKit = reorganize(fullMethylKitObj, 
#                          sample.ids=myMetaData$ID,
#                          treatment=myMetaData$trtG1G2_NUM)
# 
# ### Calculate DMS accounting for covariates: family
# cov = data.frame(Family = myMetaData$Family)
# myDiffMeth=calculateDiffMeth(myMethylKit, covariates = cov, mc.cores = 4)
# 
# # We select the bases that have q-value<0.01 and percent methylation difference larger than 15%.
# # NB: arg type="hyper" or type="hypo" gives hyper-methylated or hypo-methylated regions/bases.
# myDiff2_15p = getMethylDiff(myDiffMeth,difference=15,qvalue=0.01)
# 
# myDiff2_15p # 66 positions
# 
# ##########################################################


