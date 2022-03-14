## Differential methylation analyses
## A. Balard
## February 2022

#################### Data load & preparation ####################
source("librariesLoading.R")
## load custom functions
source("customRfunctions.R")
## Load samples metadata
source("R02.1_loadMetadata.R")
## define in which machine we're working (apocrita or mythinkpad)
##machine="apocrita"
machine="mythinkpad"
## Load methylation data
source("R02.2_loadMethyldata.R")

###########################################
## Source the previously calculated DMS/DMR
## Parents (Family as covariates)
### DM from CpG positions shared by half the fish per trt
DMS15pc_G1_half <- readRDS("../../data/DiffMeth/DMS15pc_G1_half.RDS"); nrow(DMS15pc_G1_half) # 5074
DMR15pc_G1_half <- readRDS("../../data/DiffMeth/DMR15pc_G1_half.RDS"); nrow(DMR15pc_G1_half) # 23
### DM from CpG positions shared by all the fish
DMS15pc_G1_ALL <- readRDS("../../data/DiffMeth/DMS15pc_G1_ALL.RDS"); nrow(DMS15pc_G1_ALL) # 127
# DMR15pc_G1_ALL returned 0 DMR

## Offspring (Family & Sex as covariates)
## Control G1 - G2(trt vs control)
### DM from CpG positions shared by half the fish per trt
DMS15pc_G2_controlG1_half <- readRDS("../../data/DiffMeth/DMS15pc_G2_controlG1_half.RDS"); nrow(DMS15pc_G2_controlG1_half) # 1430
DMR15pc_G2_controlG1_half <- readRDS("../../data/DiffMeth/DMR15pc_G2_controlG1_half.RDS"); nrow(DMR15pc_G2_controlG1_half) # 6
### DM from CpG positions shared by all the fish
DMS15pc_G2_controlG1_ALL <- readRDS("../../data/DiffMeth/DMS15pc_G2_controlG1_ALL.RDS"); nrow(DMS15pc_G2_controlG1_ALL) # 38
# DMR15pc_G2_controlG1_ALL returned 0 DMR

## Infected G1 - G2(trt vs control)
### DM from CpG positions shared by half the fish per trt
DMS15pc_G2_infectedG1_half <- readRDS("../../data/DiffMeth/DMS15pc_G2_infectedG1_half.RDS"); nrow(DMS15pc_G2_infectedG1_half) # 777
DMR15pc_G2_infectedG1_half <- readRDS("../../data/DiffMeth/DMR15pc_G2_infectedG1_half.RDS"); nrow(DMR15pc_G2_infectedG1_half) # 8
### DM from CpG positions shared by all the fish
DMS15pc_G2_infectedG1_ALL <- readRDS("../../data/DiffMeth/DMS15pc_G2_infectedG1_ALL.RDS"); nrow(DMS15pc_G2_infectedG1_ALL) # 22
DMR15pc_G2_infectedG1_ALL <- readRDS("../../data/DiffMeth/DMR15pc_G2_infectedG1_ALL.RDS"); nrow(DMR15pc_G2_infectedG1_ALL) # 1
###########################################

###########################
## Function to get DMS info
myDMSinfo <- function(DMSobject, fromUniteCov){
  DMS = paste(DMSobject$chr, DMSobject$start, DMSobject$end)
  meth.diff = DMSobject$meth.diff
  direction = ifelse(DMSobject$meth.diff > 0, "hyper", "hypo")
  percentDMS = length(DMS)/nrow(fromUniteCov)*100
  return(list(DMS = DMS, meth.diff = meth.diff, direction = direction, percentDMS = percentDMS))
}

## Run the function
DMS_info_G1 <- myDMSinfo(DMS15pc_G1_half, uniteCov6_G1_woSexAndUnknowChrOVERLAP)
DMS_info_G2_G1c_final <- myDMSinfo(DMS15pc_G2_controlG1_half, uniteCov14_G2_woSexAndUnknowChrOVERLAP)
DMS_info_G2_G1i_final <- myDMSinfo(DMS15pc_G2_infectedG1_half,uniteCov14_G2_woSexAndUnknowChrOVERLAP)

## NB Kostas' results: "We found a total of 1,973 CpG sites out of 1,172,887 CpGs (0.17%)
# across the genome that showed at least 15% differential fractional methylation 
# (differentially methylated site [DMS]; q < 0.01) between infected and uninfected fish"

## Here: number of CpG sites
nrow(uniteCov14_G2_woSexAndUnknowChrOVERLAP) # 1,001,880

## Parents comparison:
length(DMS_info_G1$DMS)# 5024 DMS
DMS_info_G1$percentDMS # 0.5% of the CpGs are DMS

## Offspring from control parents comparison:
length(DMS_info_G2_G1c_final$DMS) # 1430 DMS
DMS_info_G2_G1c_final$percentDMS # 0.14% of the CpGs are DMS

## Offspring from infected parents comparison:
length(DMS_info_G2_G1i_final$DMS) # 777 DMS
DMS_info_G2_G1i_final$percentDMS # 0.08% of the CpGs are DMS

###############################
## Just considering intersecting CpGs:
myVennFUN <- function(A, B, C, catnames, myCols = c("grey","green","red")){
  futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")#to rm log files
  Venn <- venn.diagram(
    x = list(A, B, C), category.names = catnames, filename = NULL,
    margin = 0, lwd = 2, lty = 'blank', fill = myCols,
    cex = .4, fontface = "bold",fontfamily = "sans", print.mode=c("raw","percent"),
    cat.cex = 0.4, cat.fontface = "bold", cat.default.pos = "outer",
    cat.col = myCols, cat.pos = c(-27, 27, 135), cat.dist = c(0.055, 0.055, 0.055),
    cat.fontfamily = "sans", rotation = 1
  )
  return(Venn)
}

## output Venn diagrams
allVenn <- myVennFUN(A = DMS_info_G1$DMS, B = DMS_info_G2_G1c_final$DMS, C = DMS_info_G2_G1i_final$DMS,
                     catnames = c("DMS G1" , "DMS G2-c", "DMS G2-i"))
hypoVenn <- myVennFUN(A = DMS_info_G1$DMS[DMS_info_G1$direction %in% "hypo"],
                      B = DMS_info_G2_G1c_final$DMS[DMS_info_G2_G1c_final$direction %in% "hypo"],
                      C = DMS_info_G2_G1i_final$DMS[DMS_info_G2_G1i_final$direction %in% "hypo"],
                      catnames = c("DMS G1\nhypo" , "DMS G2-c\nhypo", "DMS G2-i\nhypo"))
hyperVenn <- myVennFUN(A = DMS_info_G1$DMS[DMS_info_G1$direction %in% "hyper"],
                       B = DMS_info_G2_G1c_final$DMS[DMS_info_G2_G1c_final$direction %in% "hyper"],
                       C = DMS_info_G2_G1i_final$DMS[DMS_info_G2_G1i_final$direction %in% "hyper"],
                       catnames = c("DMS G1\nhyper" , "DMS G2-c\nhyper", "DMS G2-i\nhyper"))

# Output the diagrams
png(file="Rfigures/VennDMSinhalffish_intersectingCpGs_hyper_hypo.png", width = 5, height = 5.5, units = 'in', res = 300)
pushViewport(plotViewport(layout=grid.layout(2, 2)))
pushViewport(plotViewport(layout.pos.col=1, layout.pos.row=1))
grid.draw(allVenn)
popViewport()
pushViewport(plotViewport(layout.pos.col=1, layout.pos.row=2))
grid.draw(hypoVenn)
popViewport()
pushViewport(plotViewport(layout.pos.col=2, layout.pos.row=2))
grid.draw(hyperVenn)
dev.off()

rm(allVenn, hypoVenn, hyperVenn)

######################
## Features Annotation (use package genomation v1.24.0)

## load genome annotation
gene.obj=readTranscriptFeatures("../../gitignore/bigdata/Gy_allnoM_rd3.maker_apocrita.noseq_corrected.bed12", 
                                remove.unusual = FALSE)

## NB Promoters are defined by options at genomation::readTranscriptFeatures function. The default option is to take -1000,+1000bp around the TSS and you can change that. S

## annotate differentially methylated CpGs with promoter/exon/intron using annotation data
## Kostas MBE: The DMSs and regions were predominately
# found in intergenic regions (47.74% and 48.94%, respectively),
# with introns (26.19% and 23.09), exons (15.07% and 13.98%), and promoters (11% and 13.98%) showing lower proportions

par(mfrow=c(1,3))
par(mar = c(.1,0.1,5,0.1)) # Set the margin on all sides to 2
## Parents comparison:
diffAnn_PAR = annotateWithGeneParts(as(DMS15pc_G1_half,"GRanges"),gene.obj)
diffAnn_PAR
plotTargetAnnotation(diffAnn_PAR,precedence=TRUE, main="DMS G1", 
                     cex.legend = 1, border="white")
## Offspring from control parents comparison:
diffAnn_G2_controlG1 = annotateWithGeneParts(as(DMS15pc_G2_controlG1_half,"GRanges"),gene.obj)
diffAnn_G2_controlG1
plotTargetAnnotation(diffAnn_G2_controlG1,precedence=TRUE, main="DMS G2-G1c", 
                     cex.legend = 1, border="white")
## Offspring from infected parents comparison:
diffAnn_G2_infectedG1 = annotateWithGeneParts(as(DMS15pc_G2_infectedG1_half,"GRanges"),gene.obj)
diffAnn_G2_infectedG1
plotTargetAnnotation(diffAnn_G2_infectedG1,precedence=TRUE, main="DMS G2-G1i", 
                     cex.legend = 1, border="white")
par(mfrow=c(1,1)) 

## Get distance to TSS:
head(getAssociationWithTSS(diffAnn_PAR))
head(diffAnn_PAR@members)

test <- head(gene.obj)
test$TSSes

##########################
## Separate hyper and hypo
runHyperHypoAnnot <- function(){
  par(mfrow=c(2,3))
  par(mar = c(.1,0.1,5,0.1)) # Set the margin on all sides to 2
  ####### HYPO
  ## Parents comparison:
  A = annotateWithGeneParts(
    as(DMS15pc_G1_half[DMS_info_G1$direction %in% "hypo",],"GRanges"),gene.obj)
  plotTargetAnnotation(A,precedence=TRUE, main="DMS G1\nhypo", 
                       cex.legend = .4, border="white")
  ## Offspring from control parents comparison:
  B = annotateWithGeneParts(
    as(DMS15pc_G2_controlG1_half[DMS_info_G2_G1c_final$direction %in% "hypo",],"GRanges"),gene.obj)
  plotTargetAnnotation(B,precedence=TRUE, main="DMS G2-G1c\nhypo", 
                       cex.legend = .4, border="white")
  ## Offspring from infected parents comparison:
  C = annotateWithGeneParts(
    as(DMS15pc_G2_infectedG1_half[DMS_info_G2_G1i_final$direction %in% "hypo",],"GRanges"),gene.obj)
  plotTargetAnnotation(C,precedence=TRUE, main="DMS G2-G1i\nhypo", 
                       cex.legend = .4, border="white")
  ####### HYPER
  ## Parents comparison:
  D = annotateWithGeneParts(
    as(DMS15pc_G1_half[DMS_info_G1$direction %in% "hyper",],"GRanges"),gene.obj)
  plotTargetAnnotation(D,precedence=TRUE, main="DMS G1\nhyper", 
                       cex.legend = .4, border="white")
  ## Offspring from control parents comparison:
  E = annotateWithGeneParts(
    as(DMS15pc_G2_controlG1_half[DMS_info_G2_G1c_final$direction %in% "hyper",],"GRanges"),gene.obj)
  plotTargetAnnotation(E,precedence=TRUE, main="DMS G2-G1c\nhyper", 
                       cex.legend = .4, border="white")
  ## Offspring from infected parents comparison:
  f = annotateWithGeneParts(
    as(DMS15pc_G2_infectedG1_half[DMS_info_G2_G1i_final$direction %in% "hyper",],"GRanges"),gene.obj)
  plotTargetAnnotation(f,precedence=TRUE, main="DMS G2-G1i\nhyper", 
                       cex.legend = .4, border="white")
  par(mfrow=c(1,1))
  return(list(G1hypo=A, G2G1chypo=B, G2G1ihypo=C, G1hyper=D, G2G1chyper=E, G2G1ihyper=f))
}

myannot=runHyperHypoAnnot()

############################################################
## Venn diagram of overlapping features by their annotation:
table(rowSums(as.data.frame(myannot$G1hypo@members))) # NB: some positions are labelled with several features!
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

DMS15pc_G1_half_HYPO = myAnnotateDMS(DMS15pc_G1_half[DMS_info_G1$direction %in% "hypo",],
                                  as.data.frame(myannot$G1hypo@members))
DMS15pc_G1_half_HYPER = myAnnotateDMS(DMS15pc_G1_half[DMS_info_G1$direction %in% "hyper",],
                                   as.data.frame(myannot$G1hyper@members))

DMS15pc_G2_controlG1_half_HYPO = myAnnotateDMS(DMS15pc_G2_controlG1_half[DMS_info_G2_G1c_final$direction %in% "hypo",],
                                      as.data.frame(myannot$G2G1chypo@members))
DMS15pc_G2_controlG1_half_HYPER = myAnnotateDMS(DMS15pc_G2_controlG1_half[DMS_info_G2_G1c_final$direction %in% "hyper",],
                                       as.data.frame(myannot$G2G1chyper@members))

DMS15pc_G2_infectedG1_half_HYPO = myAnnotateDMS(DMS15pc_G2_infectedG1_half[DMS_info_G2_G1i_final$direction %in% "hypo",],
                                      as.data.frame(myannot$G2G1ihypo@members))
DMS15pc_G2_infectedG1_half_HYPER = myAnnotateDMS(DMS15pc_G2_infectedG1_half[DMS_info_G2_G1i_final$direction %in% "hyper",],
                                       as.data.frame(myannot$G2G1ihyper@members))

## Make Venn diagram for each feature
getFeatureDFHYPO <- function(myfeat){
  a = DMS15pc_G1_half_HYPO$pos[DMS15pc_G1_half_HYPO$feature %in% myfeat]
  b = DMS15pc_G2_controlG1_half_HYPO$pos[DMS15pc_G2_controlG1_half_HYPO$feature %in% myfeat]
  c = DMS15pc_G2_infectedG1_half_HYPO$pos[DMS15pc_G2_infectedG1_half_HYPO$feature %in% myfeat]
  return(list(a=a,b=b,c=c))
}

getFeatureDFHYPER <- function(myfeat){
  a = DMS15pc_G1_half_HYPER$pos[DMS15pc_G1_half_HYPER$feature %in% myfeat]
  b = DMS15pc_G2_controlG1_half_HYPER$pos[DMS15pc_G2_controlG1_half_HYPER$feature %in% myfeat]
  c = DMS15pc_G2_infectedG1_half_HYPER$pos[DMS15pc_G2_infectedG1_half_HYPER$feature %in% myfeat]
  return(list(a=a,b=b,c=c))
}

getVenn <- function(feat, direction){
  if (direction == "hypo"){
    myVennFUN(getFeatureDFHYPO(feat)[["a"]],
              getFeatureDFHYPO(feat)[["b"]],
              getFeatureDFHYPO(feat)[["c"]],
              catnames = c(paste0("DMS G1\nhypo\n", feat), 
                           paste0("DMS G2-c\nhypo\n", feat),
                           paste0("DMS G2-i\nhypo\n", feat)))
  }else if (direction == "hyper"){
    myVennFUN(getFeatureDFHYPER(feat)[["a"]],
              getFeatureDFHYPER(feat)[["b"]],
              getFeatureDFHYPER(feat)[["c"]],
              catnames = c(paste0("DMS G1\nhyper\n", feat), 
                           paste0("DMS G2-c\nhyper\n", feat),
                           paste0("DMS G2-i\nhyper\n", feat)))
  }
}

## output Venn diagrams
png(file="Rfigures/VennDMSinhalffish_intersectingCpGs_byfeatures_HYPO.png", width = 5, height = 5.5, units = 'in', res = 300)
#height = 500, width = 500, compression = "lzw")
pushViewport(plotViewport(layout=grid.layout(2, 2)))
pushViewport(plotViewport(layout.pos.col=1, layout.pos.row=1))
grid.draw(getVenn("promoter", "hypo"))
popViewport()
pushViewport(plotViewport(layout.pos.col=2, layout.pos.row=1))
grid.draw(getVenn("exon", "hypo"))
popViewport()
pushViewport(plotViewport(layout.pos.col=1, layout.pos.row=2))
grid.draw(grid.draw(getVenn("intron", "hypo")))
popViewport()
pushViewport(plotViewport(layout.pos.col=2, layout.pos.row=2))
grid.draw(grid.draw(getVenn("intergenic", "hypo"))) 
dev.off()

png(file="Rfigures/VennDMSinhalffish_intersectingCpGs_byfeatures_HYPER.png", width = 5, height = 5.5, units = 'in', res = 300)
#height = 500, width = 500, compression = "lzw")
pushViewport(plotViewport(layout=grid.layout(2, 2)))
pushViewport(plotViewport(layout.pos.col=1, layout.pos.row=1))
grid.draw(getVenn("promoter", "hyper"))
popViewport()
pushViewport(plotViewport(layout.pos.col=2, layout.pos.row=1))
grid.draw(getVenn("exon", "hyper"))
popViewport()
pushViewport(plotViewport(layout.pos.col=1, layout.pos.row=2))
grid.draw(grid.draw(getVenn("intron", "hyper")))
popViewport()
pushViewport(plotViewport(layout.pos.col=2, layout.pos.row=2))
grid.draw(grid.draw(getVenn("intergenic", "hyper"))) 
dev.off()

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
makeManhattanPlots(DMSfile = DMS15pc_G1_half, annotFile = annot_PAR, GYgynogff = GYgynogff, 
                   mycols = c("red", "grey", "black", "green"), mytitle = "Manhattan plot of G1 DMS")

## G2-G1c trt-ctrl
# load annotation
annot_G2_G1c <- as.data.frame(diffAnn_G2_controlG1@members)
makeManhattanPlots(DMSfile = DMS15pc_G2_controlG1_half, annotFile = annot_G2_G1c, GYgynogff = GYgynogff, 
                   mycols = c("red", "grey", "black", "green"), mytitle = "Manhattan plot of G2-G1c DMS")

## G2-G1i trt-ctrl
# load annotation
annot_G2_G1i <- as.data.frame(diffAnn_G2_infectedG1@members)
makeManhattanPlots(DMSfile = DMS15pc_G2_infectedG1_half, annotFile = annot_G2_G1i, GYgynogff = GYgynogff, 
                   mycols = c("red", "grey", "black", "green"), mytitle = "Manhattan plot of G2-G1i DMS")

## Outliers in Manhattan plot: 15% diff + 2SD
outliers_G1_final <- which(abs(DMS15pc_G1_half$meth.diff) > 15 + 2*sd(abs(DMS15pc_G1_half$meth.diff)))
outliers_annot_G1 <- as.data.frame(diffAnn_PAR@members)[outliers_G1_final,]
makeManhattanPlots(DMSfile = DMS15pc_G1_half[outliers_G1_final, ],
                   annotFile = outliers_annot_G1, GYgynogff = GYgynogff, 
                   mycols = c("red", "grey", "black", "green"), mytitle = "Manhattan plot of G1 DMS")

outliers_G2_G1c_final <- which(abs(DMS15pc_G2_controlG1_half$meth.diff) > 15 + 2*sd(abs(DMS15pc_G2_controlG1_half$meth.diff)))
outliers_annot_G2_G1c <- as.data.frame(diffAnn_G2_controlG1@members)[outliers_G2_G1c_final,]
makeManhattanPlots(DMSfile = DMS15pc_G2_controlG1_half[outliers_G2_G1c_final, ],
                   annotFile = outliers_annot_G2_G1c, GYgynogff = GYgynogff, 
                   mycols = c("red", "grey", "black", "green"), mytitle = "Manhattan plot of G2-G1c DMS")

outliers_G2_G1i_final <- which(abs(DMS15pc_G2_infectedG1_half$meth.diff) > 15 + 2*sd(abs(DMS15pc_G2_infectedG1_half$meth.diff)))
outliers_annot_G2_G1i <- as.data.frame(diffAnn_G2_infectedG1@members)[outliers_G2_G1i_final,]
makeManhattanPlots(DMSfile = DMS15pc_G2_infectedG1_half[outliers_G2_G1i_final, ],
                   annotFile = outliers_annot_G2_G1i, GYgynogff = GYgynogff, 
                   mycols = c("red", "grey", "black", "green"), mytitle = "Manhattan plot of G2-G1i DMS")

