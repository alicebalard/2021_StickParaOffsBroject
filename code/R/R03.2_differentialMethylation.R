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
## machine="apocrita"
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

## Calculate DMS accounting for covariates: family and sex (new 20/02/22!)
getDMS <- function(myuniteCov, myMetadata){
  if (length(table(myMetadata$Sex)) == 1){
    cov = data.frame(Family = myMetadata$Family)
  } else if (length(table(myMetadata$Sex)) == 2){
    cov = data.frame(Family = myMetadata$Family, Sex = myMetadata$Sex)
  } 
  myDiffMeth=calculateDiffMeth(myuniteCov, covariates = cov, mc.cores = 10)
  ## We select the bases that have q-value<0.01 and percent methylation difference larger than 15%.
  ## NB: arg type="hyper" or type="hypo" gives hyper-methylated or hypo-methylated regions/bases.
  myDMS_15pc = getMethylDiff(myDiffMeth, difference=15, qvalue=0.01)
  return(myDMS_15pc)
}

## Correspondance trtG1G2 and numerical values used by MethylKit
table(fullMetadata$trtG1G2_NUM, fullMetadata$trtG1G2)
# Control Exposed NE_control NE_exposed E_control E_exposed
# 1      12       0          0          0         0         0
# 2       0       0          0          0        28         0
# 3       0       0          0          0         0        28
# 4       0      12          0          0         0         0
# 5       0       0         28          0         0         0
# 6       0       0          0         27         0         0

## Comparison 1: BASELINE -> Parents (control vs infected) 
# DMS15pc_PAR_half <- getDMS(uniteCov6_G1_woSexAndUnknowChr, fullMetadata_PAR)
# saveRDS(DMS15pc_PAR_half, file = "../../data/DMS15pc_PAR_half.RDS")

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
# saveRDS(DMS15pc_G2_controlG1_half, file = "../../data/DMS15pc_G2_controlG1_half.RDS")
# 
# DMS15pc_G2_infectedG1_half <- getDMS(myuniteCov = reorganize(methylObj = uniteCov14_G2_woSexAndUnknowChr,
#                                                              treatment = fullMetadata_OFFS$trtG1G2_NUM[
#                                                                                                fullMetadata_OFFS$trtG1G2_NUM %in% c(2,3)], 
#                                                              sample.ids = fullMetadata_OFFS$ID[
#                                                                                                 fullMetadata_OFFS$trtG1G2_NUM %in% c(2,3)]), 
#                                      myMetadata = fullMetadata_OFFS[fullMetadata_OFFS$trtG1G2_NUM %in% c(2,3),])
#  
# saveRDS(DMS15pc_G2_infectedG1_half, file = "../../data/DMS15pc_G2_infectedG1_half.RDS")

## stop here:
# stop("We stop here for now") # to run getDMS on Apocrita cause it's LONG

## To do on bash: rename with DATES 
## for f in DMS15pc*; do mv "$f" "$(echo "$f" | sed s/.RDS/_21feb22.RDS/)"; done

### UPLOAD methylDiff objects:
# Parents (Family as covariates)
DMS15pc_PAR_half <- readRDS("../../data/DMS15pc_PAR_half_21feb22.RDS")
## Offspring (Family & Sex as covariates)
## Control G1 - G2(trt vs control)
DMS15pc_G2_controlG1_half <- readRDS("../../data/DMS15pc_G2_controlG1_half_21feb22.RDS")
## Infected G1 - G2(trt vs control)
DMS15pc_G2_infectedG1_half <- readRDS("../../data/DMS15pc_G2_infectedG1_half_21feb22.RDS")

##########################################################
##### Get DMS after removal of CpG NOT covered in all 3 datasets (G1 & G2c&i)
CpG_G1df <- paste(uniteCov6_G1_woSexAndUnknowChr$chr, uniteCov6_G1_woSexAndUnknowChr$start, uniteCov6_G1_woSexAndUnknowChr$end)
CpG_G2df <- paste(uniteCov14_G2_woSexAndUnknowChr$chr, uniteCov14_G2_woSexAndUnknowChr$start, uniteCov14_G2_woSexAndUnknowChr$end)

## Full CpG covered
coveredCpGbothdf <- intersect(CpG_G1df, CpG_G2df)
length(coveredCpGbothdf) # 1,001,880

## Keep only these positions in the 3 DMS datasets:
DMS_G1_final <- DMS15pc_PAR_half[paste(DMS15pc_PAR_half$chr, DMS15pc_PAR_half$start, DMS15pc_PAR_half$end) %in% coveredCpGbothdf,]
DMS_G2_G1c_final <- DMS15pc_G2_controlG1_half[paste(DMS15pc_G2_controlG1_half$chr, DMS15pc_G2_controlG1_half$start, DMS15pc_G2_controlG1_half$end) %in% coveredCpGbothdf,]
DMS_G2_G1i_final <- DMS15pc_G2_infectedG1_half[paste(DMS15pc_G2_infectedG1_half$chr, DMS15pc_G2_infectedG1_half$start, DMS15pc_G2_infectedG1_half$end) %in% coveredCpGbothdf,]

rm(DMS15pc_PAR_half, DMS15pc_G2_controlG1_half, DMS15pc_G2_infectedG1_half)

## Function to get DMS info
myDMSinfo <- function(DMSobject){
  DMS = paste(DMSobject$chr, DMSobject$start, DMSobject$end)
  meth.diff = DMSobject$meth.diff
  direction = ifelse(DMSobject$meth.diff > 0, "hyper", "hypo")
  percentDMS = length(DMS)/length(coveredCpGbothdf)*100
  return(list(DMS = DMS, meth.diff = meth.diff, direction = direction, percentDMS = percentDMS))
}

## Run the function (takes a few minutes)
DMS_info_G1 <- myDMSinfo(DMS_G1_final)
DMS_info_G2_G1c_final <- myDMSinfo(DMS_G2_G1c_final)
DMS_info_G2_G1i_final <- myDMSinfo(DMS_G2_G1i_final)

## NB Kostas' results: "We found a total of 1,973 CpG sites out of 1,172,887 CpGs (0.17%) across the genome that showed at
# least 15% differential fractional methylation (differentially methylated site [DMS]; q < 0.01) between infected and uninfected fish"

## Here: number of CpG sites
length(coveredCpGbothdf) # 1,001,880

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

# Outpur the diagrams
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

######################
## Features Annotation (use package genomation v1.24.0)

## load genome annotation
gene.obj=readTranscriptFeatures("../../gitignore/bigdata/Gy_allnoM_rd3.maker_apocrita.noseq_corrected.bed12", 
                                remove.unusual = FALSE)

## annotate differentially methylated CpGs with promoter/exon/intron using annotation data
## Kostas MBE: The DMSs and regions were predominately
# found in intergenic regions (47.74% and 48.94%, respectively),
# with introns (26.19% and 23.09), exons (15.07% and 13.98%), and promoters (11% and 13.98%) showing lower proportions

par(mfrow=c(1,3))
par(mar = c(.1,0.1,5,0.1)) # Set the margin on all sides to 2
## Parents comparison:
diffAnn_PAR = annotateWithGeneParts(as(DMS_G1_final,"GRanges"),gene.obj)
diffAnn_PAR
plotTargetAnnotation(diffAnn_PAR,precedence=TRUE, main="DMS G1", 
                     cex.legend = 1, border="white")
## Offspring from control parents comparison:
diffAnn_G2_controlG1 = annotateWithGeneParts(as(DMS_G2_G1c_final,"GRanges"),gene.obj)
diffAnn_G2_controlG1
plotTargetAnnotation(diffAnn_G2_controlG1,precedence=TRUE, main="DMS G2-G1c", 
                     cex.legend = 1, border="white")
## Offspring from infected parents comparison:
diffAnn_G2_infectedG1 = annotateWithGeneParts(as(DMS_G2_G1i_final,"GRanges"),gene.obj)
diffAnn_G2_infectedG1
plotTargetAnnotation(diffAnn_G2_infectedG1,precedence=TRUE, main="DMS G2-G1i", 
                     cex.legend = 1, border="white")
par(mfrow=c(1,1)) 

##########################
## Separate hyper and hypo

runHyperHypoAnnot <- function(){
  par(mfrow=c(2,3))
  par(mar = c(.1,0.1,5,0.1)) # Set the margin on all sides to 2
  ####### HYPO
  ## Parents comparison:
  A = annotateWithGeneParts(
    as(DMS_G1_final[DMS_info_G1$direction %in% "hypo",],"GRanges"),gene.obj)
  plotTargetAnnotation(A,precedence=TRUE, main="DMS G1\nhypo", 
                       cex.legend = .4, border="white")
  ## Offspring from control parents comparison:
  B = annotateWithGeneParts(
    as(DMS_G2_G1c_final[DMS_info_G2_G1c_final$direction %in% "hypo",],"GRanges"),gene.obj)
  plotTargetAnnotation(B,precedence=TRUE, main="DMS G2-G1c\nhypo", 
                       cex.legend = .4, border="white")
  ## Offspring from infected parents comparison:
  C = annotateWithGeneParts(
    as(DMS_G2_G1i_final[DMS_info_G2_G1i_final$direction %in% "hypo",],"GRanges"),gene.obj)
  plotTargetAnnotation(C,precedence=TRUE, main="DMS G2-G1i\nhypo", 
                       cex.legend = .4, border="white")
  ####### HYPER
  ## Parents comparison:
  D = annotateWithGeneParts(
    as(DMS_G1_final[DMS_info_G1$direction %in% "hyper",],"GRanges"),gene.obj)
  plotTargetAnnotation(D,precedence=TRUE, main="DMS G1\nhyper", 
                       cex.legend = .4, border="white")
  ## Offspring from control parents comparison:
  E = annotateWithGeneParts(
    as(DMS_G2_G1c_final[DMS_info_G2_G1c_final$direction %in% "hyper",],"GRanges"),gene.obj)
  plotTargetAnnotation(E,precedence=TRUE, main="DMS G2-G1c\nhyper", 
                       cex.legend = .4, border="white")
  ## Offspring from infected parents comparison:
  f = annotateWithGeneParts(
    as(DMS_G2_G1i_final[DMS_info_G2_G1i_final$direction %in% "hyper",],"GRanges"),gene.obj)
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

DMS_G1_final_HYPO = myAnnotateDMS(DMS_G1_final[DMS_info_G1$direction %in% "hypo",],
                                  as.data.frame(myannot$G1hypo@members))
DMS_G1_final_HYPER = myAnnotateDMS(DMS_G1_final[DMS_info_G1$direction %in% "hyper",],
                                   as.data.frame(myannot$G1hyper@members))

DMS_G2_G1c_final_HYPO = myAnnotateDMS(DMS_G2_G1c_final[DMS_info_G2_G1c_final$direction %in% "hypo",],
                                      as.data.frame(myannot$G2G1chypo@members))
DMS_G2_G1c_final_HYPER = myAnnotateDMS(DMS_G2_G1c_final[DMS_info_G2_G1c_final$direction %in% "hyper",],
                                       as.data.frame(myannot$G2G1chyper@members))

DMS_G2_G1i_final_HYPO = myAnnotateDMS(DMS_G2_G1i_final[DMS_info_G2_G1i_final$direction %in% "hypo",],
                                      as.data.frame(myannot$G2G1ihypo@members))
DMS_G2_G1i_final_HYPER = myAnnotateDMS(DMS_G2_G1i_final[DMS_info_G2_G1i_final$direction %in% "hyper",],
                                       as.data.frame(myannot$G2G1ihyper@members))

## Make Venn diagram for each feature
getFeatureDFHYPO <- function(myfeat){
  a = DMS_G1_final_HYPO$pos[DMS_G1_final_HYPO$feature %in% myfeat]
  b = DMS_G2_G1c_final_HYPO$pos[DMS_G2_G1c_final_HYPO$feature %in% myfeat]
  c = DMS_G2_G1i_final_HYPO$pos[DMS_G2_G1i_final_HYPO$feature %in% myfeat]
  return(list(a=a,b=b,c=c))
}

getFeatureDFHYPER <- function(myfeat){
  a = DMS_G1_final_HYPER$pos[DMS_G1_final_HYPER$feature %in% myfeat]
  b = DMS_G2_G1c_final_HYPER$pos[DMS_G2_G1c_final_HYPER$feature %in% myfeat]
  c = DMS_G2_G1i_final_HYPER$pos[DMS_G2_G1i_final_HYPER$feature %in% myfeat]
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
makeManhattanPlots(DMSfile = DMS_G1_final, annotFile = annot_PAR, GYgynogff = GYgynogff, 
                   mycols = c("red", "grey", "black", "green"), mytitle = "Manhattan plot of G1 DMS")

## G2-G1c trt-ctrl
# load annotation
annot_G2_G1c <- as.data.frame(diffAnn_G2_controlG1@members)
makeManhattanPlots(DMSfile = DMS_G2_G1c_final, annotFile = annot_G2_G1c, GYgynogff = GYgynogff, 
                   mycols = c("red", "grey", "black", "green"), mytitle = "Manhattan plot of G2-G1c DMS")

## G2-G1i trt-ctrl
# load annotation
annot_G2_G1i <- as.data.frame(diffAnn_G2_infectedG1@members)
makeManhattanPlots(DMSfile = DMS_G2_G1i_final, annotFile = annot_G2_G1i, GYgynogff = GYgynogff, 
                   mycols = c("red", "grey", "black", "green"), mytitle = "Manhattan plot of G2-G1i DMS")

## Outliers in Manhattan plot: 15% diff + 2SD
outliers_G1_final <- which(abs(DMS_G1_final$meth.diff) > 15 + 2*sd(abs(DMS_G1_final$meth.diff)))
outliers_annot_G1 <- as.data.frame(diffAnn_PAR@members)[outliers_G1_final,]
makeManhattanPlots(DMSfile = DMS_G1_final[outliers_G1_final, ],
                   annotFile = outliers_annot_G1, GYgynogff = GYgynogff, 
                   mycols = c("red", "grey", "black", "green"), mytitle = "Manhattan plot of G1 DMS")

outliers_G2_G1c_final <- which(abs(DMS_G2_G1c_final$meth.diff) > 15 + 2*sd(abs(DMS_G2_G1c_final$meth.diff)))
outliers_annot_G2_G1c <- as.data.frame(diffAnn_G2_controlG1@members)[outliers_G2_G1c_final,]
makeManhattanPlots(DMSfile = DMS_G2_G1c_final[outliers_G2_G1c_final, ],
                   annotFile = outliers_annot_G2_G1c, GYgynogff = GYgynogff, 
                   mycols = c("red", "grey", "black", "green"), mytitle = "Manhattan plot of G2-G1c DMS")

outliers_G2_G1i_final <- which(abs(DMS_G2_G1i_final$meth.diff) > 15 + 2*sd(abs(DMS_G2_G1i_final$meth.diff)))
outliers_annot_G2_G1i <- as.data.frame(diffAnn_G2_infectedG1@members)[outliers_G2_G1i_final,]
makeManhattanPlots(DMSfile = DMS_G2_G1i_final[outliers_G2_G1i_final, ],
                   annotFile = outliers_annot_G2_G1i, GYgynogff = GYgynogff, 
                   mycols = c("red", "grey", "black", "green"), mytitle = "Manhattan plot of G2-G1i DMS")


############################## Identify Genes associated with positions 

############################## + GO terms --> TBC

### TBC: DMR 
# The function below tiles the genome with windows of 1000 bp length and 1000 bp
# step-size and summarizes the methylation information on those tiles. In this case, it
# returns a methylRawList object which can be fed into unite() and calculateDiffMeth()
# functions consecutively to get differentially methylated regions:

# (check Kostas and Mel papers to see which range they chose)

# tiles=tileMethylCounts(myobj,win.size=1000,step.size=1000)
# head(tiles[[1]],3)
