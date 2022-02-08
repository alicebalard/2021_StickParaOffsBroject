## Differential methylation analyses
## A. Balard
## February 2022

#################### Data load & preparation ####################
source("librariesLoading.R")
## load custom functions
source("customRfunctions.R")
## Load samples metadata
source("R01.3_prepMetadata.R")
## define in which machine we're working (apocrita or mythinkpad)
machine="apocrita"
## machine="mythinkpad"
## Load methylation data
source("R01.4_prepMethyldata.R")

###############################################
### PART 2: Differential methylation sites ####
###############################################

##########################################################
## Comparison 1: BASELINE -> Parents (control vs infected) 
## Extract DMS from parents (at least in 2 fish), annotate them, compare with Kostas results

## Calculate DMS accounting for covariates: family
rerun=FALSE
if (rerun==TRUE){
    cov = data.frame(Family = fullMetadata_PAR$Family)
    myDiffMeth=calculateDiffMeth(uniteCov2_woSexAndUnknowChr_PAR, 
                                 covariates = cov, mc.cores = 4)

                                        # We select the bases that have q-value<0.01 and percent methylation difference larger than 15%.
                                        # NB: arg type="hyper" or type="hypo" gives hyper-methylated or hypo-methylated regions/bases.
    myDiff1_15p = getMethylDiff(myDiffMeth,difference=15,qvalue=0.01)

    myDiff1_15p # 6544 positions
    saveRDS(myDiff1_15p, file = "../../data/myDiff1_15p_parentalDiffMeth.RDS")
}

DMS_PAR <- readRDS("../../data/myDiff1_15p_parentalDiffMeth.RDS")

#############
## annotation

## load genome annotation
gene.obj=readTranscriptFeatures("../../gitignore/bigdata/Gy_allnoM_rd3.maker_apocrita.noseq_corrected.bed12", remove.unusual = FALSE)
 
# annotate differentially methylated CpGs with promoter/exon/intron using annotation data
diffAnn_PAR=annotateWithGeneParts(as(DMS_PAR,"GRanges"),gene.obj)
diffAnn_PAR
plotTargetAnnotation(diffAnn_PAR,precedence=TRUE,
                     main="differential methylation annotation")
## Kostas MBE: The DMSs and regions were predominately
#found in intergenic regions (47.74% and 48.94%, respecti vely),
#with introns (26.19% and 23.09), exons (15.07% and 13.98%),and promoters (11% and 13.98%) showing lower proportions
 
###########################
## Manhattan plot of DMS ##
###########################

# load file with your DMS
DMS_PAR <- readRDS("../../data/myDiff1_15p_parentalDiffMeth.RDS")
# load annotation
annot_PAR <- as.data.frame(diffAnn_PAR@members)
## Load file containing length of each gynogen chromosomes 
## grep "contig" gitignore/bigdata/Gy_allnoM_rd3.maker_apocrita.noseq_corrected.gff | awk '{print $1, $5}' > data/Gy_allnoM_rd3.maker_apocrita.noseq_corrected_chromoAndLength.txt
GYgynogff = read.table("../../data/Gy_allnoM_rd3.maker_apocrita.noseq_corrected_chromoAndLength.txt")
names(GYgynogff) = c("chrom","length")

makeManhattanPlots(DMSfile = DMS_PAR, annotFile = annot_PAR, GYgynogff = GYgynogff)


# GO terms --> TBC
# Networks --> TBC


                                        # ### YOU'RE HERE ;)
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


