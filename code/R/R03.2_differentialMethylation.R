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

## Comparison 1: BASELINE -> Parents (control vs infected) 
# DMS15pc_PAR_half <- getDMS(uniteCov6_G1_woSexAndUnknowChr, fullMetadata_PAR_half)
# saveRDS(DMS15pc_PAR_half, file = "../../data/DMS15pc_PAR_half.RDS")
DMS15pc_PAR_half <- readRDS("../../data/DMS15pc_PAR_half.RDS")

## Comparison 2: Offspring
## 2.1. Should be like baseline -> G2 from G1 control (control vs infected) 
## 2.2. Should be DIFFERENT if there is a paternal effect -> G2 from G1 infected (control vs infected) 

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

## Infected G1 - G2(trt vs control)
DMS15pc_G2_infectedG1_half <- readRDS("../../data/DMS15pc_G2_infectedG1_half.RDS")

##########################################################
## Function to get DMS info
myDMSinfo <- function(methylObject, DMSobject){
  CpGs = paste(methylObject$chr, methylObject$start, methylObject$end)
  DMS = paste(DMSobject$chr, DMSobject$start, DMSobject$end)
  percentDMS = length(DMS)/length(CpGs)*100
  return(list(CpGs = CpGs, DMS = DMS, percentDMS = percentDMS))
}

## Run the function (takes a few minutes)
# myDMSinfo_PAR <- myDMSinfo(uniteCov6_G1_woSexAndUnknowChr, DMS15pc_PAR_half)
# myDMSinfo_G2_controlG1_half <- myDMSinfo(methylObject = reorganize(methylObj = uniteCov14_G2_woSexAndUnknowChr,
#                                                                    treatment = fullMetadata_OFFS$trtG1G2_NUM[
#                                                                      fullMetadata_OFFS$trtG1G2_NUM %in% c(5,6)],
#                                                                    sample.ids = fullMetadata_OFFS$ID[
#                                                                      fullMetadata_OFFS$trtG1G2_NUM %in% c(5,6)]),
#                                          DMSobject = DMS15pc_G2_controlG1_half)
# myDMSinfo_G2_infectedG1_half <- myDMSinfo(methylObject = reorganize(methylObj = uniteCov14_G2_woSexAndUnknowChr,
#                                                                     treatment = fullMetadata_OFFS$trtG1G2_NUM[
#                                                                       fullMetadata_OFFS$trtG1G2_NUM %in% c(2,3)],
#                                                                     sample.ids = fullMetadata_OFFS$ID[
#                                                                       fullMetadata_OFFS$trtG1G2_NUM %in% c(2,3)]),
#                                           DMSobject = DMS15pc_G2_infectedG1_half)

##### Calculate DMS after removal of CpG NOT covered in both datasets (G1 & G2)
intersectCpG <- intersect(intersect(myDMSinfo_PAR$CpGs, myDMSinfo_G2_controlG1_half$CpGs),
                          myDMSinfo_G2_infectedG1_half$CpGs)

length(intersectCpG) # 1,001,880

### Add info on DMS at intersecting CpG positions
updateIntersectFun <- function(myDMSinfo){
  myDMSinfo$DMS_intersect <- myDMSinfo$DMS[myDMSinfo$DMS %in% intersectCpG]
  myDMSinfo$percentDMS_intersect <- length(myDMSinfo$DMS_intersect) / length(intersectCpG)*100
  return(myDMSinfo)
}

##################
### Read out infos
myDMSinfo_PAR <- updateIntersectFun(myDMSinfo_PAR)
myDMSinfo_G2_controlG1_half <- updateIntersectFun(myDMSinfo_G2_controlG1_half)
myDMSinfo_G2_infectedG1_half <- updateIntersectFun(myDMSinfo_G2_infectedG1_half)

length(intersectCpG) # 1,001,880 CpG considered

## NB Kostas' results: "We found a total of 1,973 CpG sites out of 1,172,887 CpGs (0.17%) across the genome that showed at
# least 15% differential fractional methylation (differentially methylated site [DMS]; q < 0.01) between infected and uninfected fish"

## Parents comparison:
length(myDMSinfo_PAR$DMS_intersect) # 5024 DMS
myDMSinfo_PAR$percentDMS_intersect # 0.5% of the CpGs are DMS

## Offspring from control parents comparison:
length(myDMSinfo_G2_controlG1_half$DMS_intersect) # 1478 DMS
myDMSinfo_G2_controlG1_half$percentDMS_intersect # 0.15% of the CpGs are DMS

## Offspring from infected parents comparison:
length(myDMSinfo_G2_infectedG1_half$DMS_intersect) # 832 DMS
myDMSinfo_G2_infectedG1_half$percentDMS_intersect # 0.08% of the CpGs are DMS

###############################
# Plot CpG overlapping between both datasets (G1 & G2)
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")#to rm log files
myCols = c("grey","black")
venn.diagram(
  x = list(myDMSinfo_PAR$CpG, myDMSinfo_G2_controlG1_half$CpG),
  category.names = c("G1" , "G2"),
  filename = "Rfigures/VennCpGinhalffish.png", output=TRUE,
  # Output features
  imagetype="png" , height = 480 , width = 480, resolution = 300, compression = "lzw", margin = 0.05,
  # Circles
  lwd = 2, lty = 'blank', fill = myCols,
  # Numbers
  cex = .5, fontface = "bold",fontfamily = "sans", print.mode=c("raw","percent"),
  # Set names
  cat.cex = 0.4, cat.fontface = "bold", cat.default.pos = "outer",
  cat.col = myCols, cat.pos = c(-27, 27), cat.dist = c(0.015, 0.015), cat.fontfamily = "sans"
)

## Check which DMS are overlapping between the 2 half offspring and half parents datasets:
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")#to rm log files
myCols = c("grey","green","red")
venn.diagram(
  x = list(myDMSinfo_PAR$DMS, myDMSinfo_G2_controlG1_half$DMS, myDMSinfo_G2_infectedG1_half$DMS),
  category.names = c("DMS G1" , "DMS G2-c", "DMS G2-i"),
  filename = "Rfigures/VennDMSinhalffish.png", output=TRUE,
  # Output features
  imagetype="png" , height = 480 , width = 480, resolution = 300, 
  compression = "lzw", margin = 0.05,
  # Circles
  lwd = 2, lty = 'blank', fill = myCols,
  # Numbers
  cex = .4, fontface = "bold",fontfamily = "sans", print.mode=c("raw","percent"),
  # Set names
  cat.cex = 0.4, cat.fontface = "bold", cat.default.pos = "outer",
  cat.col = myCols, cat.pos = c(-27, 27, 135), cat.dist = c(0.055, 0.055, 0.055),
  cat.fontfamily = "sans", rotation = 1
)

## Just considering intersecting CpGs:
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")#to rm log files
myCols = c("grey","green","red")
venn.diagram(
  x = list(myDMSinfo_PAR$DMS_intersect, myDMSinfo_G2_controlG1_half$DMS_intersect, myDMSinfo_G2_infectedG1_half$DMS_intersect),
  category.names = c("DMS G1" , "DMS G2-c", "DMS G2-i"),
  filename = "Rfigures/VennDMSinhalffish_intersectingCpGs.png", output=TRUE,
  # Output features
  imagetype="png" , height = 480 , width = 480, resolution = 300, 
  compression = "lzw", margin = 0.05,
  # Circles
  lwd = 2, lty = 'blank', fill = myCols,
  # Numbers
  cex = .4, fontface = "bold",fontfamily = "sans", print.mode=c("raw","percent"),
  # Set names
  cat.cex = 0.4, cat.fontface = "bold", cat.default.pos = "outer",
  cat.col = myCols, cat.pos = c(-27, 27, 135), cat.dist = c(0.055, 0.055, 0.055),
  cat.fontfamily = "sans", rotation = 1
)

######################
## Features Annotation (use package genomation v1.24.0)

## load genome annotation
gene.obj=readTranscriptFeatures("../../gitignore/bigdata/Gy_allnoM_rd3.maker_apocrita.noseq_corrected.bed12", 
                                remove.unusual = FALSE)

##### Calculate DMS after removal of CpG NOT covered in both datasets (G1 & G2)
# Selects the rows corresponding to intersecting CpGs in our 3 methylDiff objects:
# DMS15pc_PAR_half
# DMS15pc_G2_controlG1_half
# DMS15pc_G2_infectedG1_half

parpos <- which(paste(DMS15pc_PAR_half$chr, DMS15pc_PAR_half$start, DMS15pc_PAR_half$end) %in% intersectCpG)
DMS15pc_PAR_half_intersect <- methylKit::select(DMS15pc_PAR_half,parpos)

G2_G1cpos <- which(paste(DMS15pc_G2_controlG1_half$chr, DMS15pc_G2_controlG1_half$start, DMS15pc_G2_controlG1_half$end) %in% intersectCpG)
DMS15pc_G2_controlG1_half_intersect <- methylKit::select(DMS15pc_G2_controlG1_half,G2_G1cpos)

G2_G1ipos <- which(paste(DMS15pc_G2_infectedG1_half$chr, DMS15pc_G2_infectedG1_half$start, DMS15pc_G2_infectedG1_half$end) %in% intersectCpG)
DMS15pc_G2_infectedG1_half_intersect <- methylKit::select(DMS15pc_G2_infectedG1_half,G2_G1ipos)

## annotate differentially methylated CpGs with promoter/exon/intron using annotation data
## Kostas MBE: The DMSs and regions were predominately
# found in intergenic regions (47.74% and 48.94%, respectively),
# with introns (26.19% and 23.09), exons (15.07% and 13.98%), and promoters (11% and 13.98%) showing lower proportions

par(mfrow=c(1,3))
par(mar = c(.1,0.1,5,0.1)) # Set the margin on all sides to 2

## Parents comparison:
diffAnn_PAR = annotateWithGeneParts(as(DMS15pc_PAR_half_intersect,"GRanges"),gene.obj)
diffAnn_PAR
plotTargetAnnotation(diffAnn_PAR,precedence=TRUE, main="DMS G1", 
                     cex.legend = 1, border="white")

## Offspring from control parents comparison:
diffAnn_G2_controlG1 = annotateWithGeneParts(as(DMS15pc_G2_controlG1_half_intersect,"GRanges"),gene.obj)
diffAnn_G2_controlG1
plotTargetAnnotation(diffAnn_G2_controlG1,precedence=TRUE, main="DMS G2-G1c", 
                     cex.legend = 1, border="white")

## Offspring from infected parents comparison:
diffAnn_G2_infectedG1 = annotateWithGeneParts(as(DMS15pc_G2_infectedG1_half_intersect,"GRanges"),gene.obj)
diffAnn_G2_infectedG1
plotTargetAnnotation(diffAnn_G2_infectedG1,precedence=TRUE, main="DMS G2-G1i", 
                     cex.legend = 1, border="white")
par(mfrow=c(1,1)) 

###########################
## Manhattan plot of DMS ##
###########################
## Load file containing length of each gynogen chromosomes 
## grep "contig" gitignore/bigdata/Gy_allnoM_rd3.maker_apocrita.noseq_corrected.gff | awk '{print $1, $5}' > data/Gy_allnoM_rd3.maker_apocrita.noseq_corrected_chromoAndLength.txt
GYgynogff = read.table("../../data/Gy_allnoM_rd3.maker_apocrita.noseq_corrected_chromoAndLength.txt")
names(GYgynogff) = c("chrom","length")

## Parents trt-ctrl
# load annotation
annot_PAR <- as.data.frame(diffAnn_PAR@members)
makeManhattanPlots(DMSfile = DMS15pc_PAR_half_intersect, annotFile = annot_PAR, GYgynogff = GYgynogff, 
                   mycols = c("red", "grey", "black", "green"), mytitle = "Manhattan plot of G1 DMS")

## G2-G1c trt-ctrl
# load annotation
annot_G2_G1c <- as.data.frame(diffAnn_G2_controlG1@members)
makeManhattanPlots(DMSfile = DMS15pc_G2_controlG1_half_intersect, annotFile = annot_G2_G1c, GYgynogff = GYgynogff, 
                   mycols = c("red", "grey", "black", "green"), mytitle = "Manhattan plot of G2-G1c DMS")

## G2-G1i trt-ctrl
# load annotation
annot_G2_G1i <- as.data.frame(diffAnn_G2_infectedG1@members)
makeManhattanPlots(DMSfile = DMS15pc_G2_infectedG1_half_intersect, annotFile = annot_G2_G1i, GYgynogff = GYgynogff, 
                   mycols = c("red", "grey", "black", "green"), mytitle = "Manhattan plot of G2-G1i DMS")

############################################################
## Venn diagram of overlapping features by their annotation:
table(rowSums(annot_PAR)) # NB: some positions are labelled with several features!
## as in MBE 2021: "giving precedence to the following order promoters, exons,
## introns, and intergenic regions when features overlapped"

myAnnotateDMS <- function(DMS, annot){
  ## sanity check
  if (nrow(DMS) != nrow(annot)){"STOP error in arguments"}
  DMS$pos <- paste(DMS$chr, DMS$start, DMS$end)
  ## NB as in MBE 2021: "giving precedence to the following order promoters, exons,
  ## introns, and intergenic regions when features overlapped"
  DMS$feature <- NA
  ## 1. promoters
  DMS$feature[which(annot$prom == 1)] = "promoter"
  ## 2. exons
  DMS$feature[which(annot$exon == 1 & annot$prom ==0)] = "exon"
  ## 3. intron
  DMS$feature[which(annot$intro == 1 & annot$exon == 0 & annot$prom ==0)] = "intron"
  ## 4. intergenic regions
  DMS$feature[which(annot$intro == 0 & annot$exon == 0 & annot$prom ==0)] = "intergenic"
  return(DMS)
}

## Annotate the DMS dataframes: 
DMS15pc_PAR_half_intersect = myAnnotateDMS(DMS15pc_PAR_half_intersect, annot_PAR)
DMS15pc_G2_controlG1_half_intersect = myAnnotateDMS(DMS15pc_G2_controlG1_half_intersect, annot_G2_G1c)
DMS15pc_G2_infectedG1_half_intersect = myAnnotateDMS(DMS15pc_G2_infectedG1_half_intersect, annot_G2_G1i)

## Make Venn diagram for each feature
myFeaturesVenn <- function(myfeat){
  a = DMS15pc_PAR_half_intersect$pos[DMS15pc_PAR_half_intersect$feature %in% myfeat]
  b = DMS15pc_G2_controlG1_half_intersect$pos[DMS15pc_G2_controlG1_half_intersect$feature %in% myfeat]
  c = DMS15pc_G2_infectedG1_half_intersect$pos[DMS15pc_G2_infectedG1_half_intersect$feature %in% myfeat]
  
  futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")#to rm log files
  myCols = c("grey","green","red")
  Venn <- venn.diagram(
    x = list(a,b,c),
    category.names = c(paste("DMS G1\n", myfeat), paste("DMS G2-G1c\n", myfeat),paste("DMS G2-G1i\n", myfeat)),
    # filename = paste0("Rfigures/VennDMSinhalffish_intersectingCpGs_", myfeat, ".png"), output=TRUE,
    # Output features
    # imagetype="png" , height = 480 , width = 480, resolution = 300, compression = "lzw", 
    margin = 0,
    filename = NULL, 
    # Circles
    lwd = 1, lty = 'blank', fill = myCols,
    # Numbers
    cex = 0.3, fontface = "bold",fontfamily = "sans", print.mode=c("raw","percent"),
    # Set names
    cat.cex = 0.4, cat.default.pos = "outer", #cat.fontface = "bold",
    cat.col = myCols, cat.pos = c(-30, 30, 180), cat.dist = c(0.15,0.15,0.15),
    cat.fontfamily = "sans", rotation = 1
  )
  return(Venn)
}

## output Venn diagrams
P <- myFeaturesVenn(myfeat = "promoter")
E <- myFeaturesVenn(myfeat = "exon")
I <- myFeaturesVenn(myfeat = "intron")
I2 <- myFeaturesVenn(myfeat = "intergenic")

# Draw the diagrams
png(file="Rfigures/VennDMSinhalffish_intersectingCpGs_byfeatures.png", width = 5, height = 5.5, units = 'in', res = 300)
#height = 500, width = 500, compression = "lzw")
pushViewport(plotViewport(layout=grid.layout(2, 2)))
pushViewport(plotViewport(layout.pos.col=1, layout.pos.row=1))
grid.draw(P)
popViewport()
pushViewport(plotViewport(layout.pos.col=2, layout.pos.row=1))
grid.draw(E)
popViewport()
pushViewport(plotViewport(layout.pos.col=1, layout.pos.row=2))
grid.draw(I)
popViewport()
pushViewport(plotViewport(layout.pos.col=2, layout.pos.row=2))
grid.draw(I2) 
dev.off()

############################## GO terms --> TBC

