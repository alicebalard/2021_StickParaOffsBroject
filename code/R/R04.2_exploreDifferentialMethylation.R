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
## Load file containing length of each gynogen chromosomes 
## grep "contig" gitignore/bigdata/Gy_allnoM_rd3.maker_apocrita.noseq_corrected.gff | awk '{print $1, $5}' > data/Gy_allnoM_rd3.maker_apocrita.noseq_corrected_chromoAndLength.txt
GYgynogff = read.table("../../data/Gy_allnoM_rd3.maker_apocrita.noseq_corrected_chromoAndLength.txt")
names(GYgynogff) = c("chrom","length")

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
length(DMS_info_G1$DMS)# 5074 DMS
DMS_info_G1$percentDMS # 0.51% of the CpGs are DMS

## Offspring from control parents comparison:
length(DMS_info_G2_G1c_final$DMS) # 1430 DMS
DMS_info_G2_G1c_final$percentDMS # 0.14% of the CpGs are DMS

## Offspring from infected parents comparison:
length(DMS_info_G2_G1i_final$DMS) # 777 DMS
DMS_info_G2_G1i_final$percentDMS # 0.08% of the CpGs are DMS

###############################################################################
## Question: how are the beta values in the 4 G2 groups for the parental DMS?##
###############################################################################

## Calculate beta values (methylation proportion per CpG site)
percMethMat_uniteCov14_G2_woSexAndUnknowChrOVERLAP = 
  methylKit::percMethylation(uniteCov14_G2_woSexAndUnknowChrOVERLAP)

## Each row is a CpG sites, let's give them a proper "pos" row name
rownames(percMethMat_uniteCov14_G2_woSexAndUnknowChrOVERLAP) <- paste(uniteCov14_G2_woSexAndUnknowChrOVERLAP$chr, 
                                                                      uniteCov14_G2_woSexAndUnknowChrOVERLAP$start, 
                                                                      uniteCov14_G2_woSexAndUnknowChrOVERLAP$end)

## Select only the positions corresponding in DMS in G1 comparison control/infected
length(DMS_info_G1$DMS)# 5074 DMS

percMethMat_halfG2_atDMS <- percMethMat_uniteCov14_G2_woSexAndUnknowChrOVERLAP[
  rownames(percMethMat_uniteCov14_G2_woSexAndUnknowChrOVERLAP) %in% DMS_info_G1$DMS, ]

nrow(percMethMat_halfG2_atDMS) # all good 

percMethMat_halfG2_atDMS <- melt(percMethMat_halfG2_atDMS)

## Extract chromosome, position, and assign correct names
percMethMat_halfG2_atDMS$Chr <- sapply(strsplit(as.character(percMethMat_halfG2_atDMS$Var1), " +"), `[`, 1)
percMethMat_halfG2_atDMS$Pos <- sapply(strsplit(as.character(percMethMat_halfG2_atDMS$Var1), " +"), `[`, 2)
names(percMethMat_halfG2_atDMS) <- c("Var1",  "ID",  "BetaValue", "Chr", "Pos")
percMethMat_halfG2_atDMS$rankpos <- 1:nrow(percMethMat_halfG2_atDMS)

## Add treatment, Sex and Family
dfTrt = data.frame(ID = fullMetadata_OFFS$SampleID, Treatment = fullMetadata_OFFS$trtG1G2,
                   Sex = fullMetadata_OFFS$Sex, Family= fullMetadata_OFFS$Family)
percMethMat_halfG2_atDMS = merge(percMethMat_halfG2_atDMS, dfTrt)
percMethMat_halfG2_atDMS$G1_trt <- sapply(strsplit(as.character(percMethMat_halfG2_atDMS$Treatment), "_"), `[`, 1)
percMethMat_halfG2_atDMS$G2_trt <- sapply(strsplit(as.character(percMethMat_halfG2_atDMS$Treatment), "_"), `[`, 2)
percMethMat_halfG2_atDMS$G1_trt[percMethMat_halfG2_atDMS$G1_trt %in% "E"] <- "exposed"
percMethMat_halfG2_atDMS$G1_trt[percMethMat_halfG2_atDMS$G1_trt %in% "NE"] <- "control"

head(percMethMat_halfG2_atDMS)

## Linear model: does the beta value at DMS depends on treatment Parent x Offspring?
modFull <- lmer(BetaValue ~ G1_trt * G2_trt + (1|Sex) + (1|Family), 
                data = percMethMat_halfG2_atDMS, REML = F) # REML =F for model comparison
mod_noG1trt <- lmer(BetaValue ~ G2_trt + (1|Sex) + (1|Family), 
                    data = percMethMat_halfG2_atDMS, REML = F)
mod_noG2trt <-lmer(BetaValue ~ G1_trt + (1|Sex) + (1|Family), 
                   data = percMethMat_halfG2_atDMS, REML = F)
mod_noInteractions <- lmer(BetaValue ~ G1_trt + G2_trt + (1|Sex) + (1|Family), 
                           data = percMethMat_halfG2_atDMS, REML = F)

lrtest(modFull, mod_noG1trt) # G1 trt is VERY VERY significant p = 4.758e-07 ***
lrtest(modFull, mod_noG2trt) # G1 trt is not significant (p<0.05)
lrtest(modFull, mod_noInteractions) # interactions are not not significant (p<0.05)

modFinal <- lmer(BetaValue ~ G1_trt + (1|Sex) + (1|Family), 
                 data = percMethMat_halfG2_atDMS)

ggpredict(modFinal, terms = c("G1_trt"))
## Higher beta values at the parental DMS for G2 from exposed G1



###############################################################################
## Beta values of offsprings at parental DMS, per trt, along the chromosomes ##
###############################################################################

## Too heavy to run on my thinkpad
# modG1G2beta <- lmer(BetaValue ~ Var1 : Treatment + (1|Sex) + (1|Family), 
#                 data = percMethMat_halfG2_atDMS, REML = T)




predPosBeta <- ggpredict(modG1G2beta, terms = c("Var1", "G1_trt", "G2_trt"))


# makeManhattanPlots <- function(DMSfile, annotFile, GYgynogff, mycols=c("grey50","grey50","darkred","darkred"), 
#                                mytitle = "Manhattan plot of DMS"){
#GA_genome.fa.sizes.txt is a file with chromosome sizes and names
genome <- GYgynogff %>%
  mutate(chrom_nr=chrom %>% deroman(),
         chrom_order=factor(chrom_nr) %>% as.numeric()) %>% 
  arrange(chrom_order) %>%
  mutate(gstart=lag(length,default=0) %>% cumsum(),
         gend=gstart+length, 
         type=LETTERS[2-(chrom_order%%2)],
         gmid=(gstart+gend)/2)

# region=as.factor(ifelse(annotFile$prom==1,"promoter",
#                         ifelse(annotFile$exon==1,"exon",
#                                ifelse(annotFile$intron==1, "intron","intergenic"))))


head(percMethMat_halfG2_atDMS)

head(melt(G2BetaAtDMS))

# DMSfile = 
head(DMS15pc_G1_half)

## Offspring methylation file:


mydata = tibble(chrom=DMSfile$chr,
                pos=DMSfile$start,
                meth.diff=DMSfile$meth.diff,
                qval=DMSfile$qvalue)#,
# region=region)

# table(DMSfile$chr)## check that chrXIX and chrUN are well removed!!

# join DMS and genomic position
data = left_join(mydata, genome2) %>% 
  mutate(gpos=pos+gstart,significance= ifelse(abs(qval>0.0125) | abs(meth.diff)<15,"not significant","significant"))

table(data$significance) # all signif

#plot only significant DMS:
ggplot()+
  geom_rect(data=genome2,aes(xmin=gstart,xmax=gend,ymin=-Inf,ymax=Inf,fill=type), alpha=.2)+
  geom_point(data=data[abs(data$meth.diff)>15 & data$significance=="significant",],
             aes(x=gpos,y=meth.diff,col=region,shape=region),fill="white", size = 2)+
  scale_color_manual(values = mycols)+
  scale_shape_manual(values=c(21,21,21,21))+
  scale_fill_manual(values=c(A=rgb(.9,.9,.9),B=NA),guide="none")+
  scale_x_continuous(breaks=genome2$gmid,labels=genome2$chrom %>% str_remove(.,"Gy_chr"),
                     position = "top",expand = c(0,0))+
  theme_minimal()+
  theme(panel.grid = element_blank(),
        axis.line=element_blank(),
        axis.title = element_blank(),
        strip.placement = "outside")+
  ggtitle(mytitle)
}












#########################
## Clustering analysis ##
#########################
## Run Adonis to this if clustering is done by treatment
x = percMethMat_uniteCov14_G2_woSexAndUnknowChrOVERLAP[
  rownames(percMethMat_uniteCov14_G2_woSexAndUnknowChrOVERLAP) %in% DMS_info_G1$DMS, ]

## Transpose
x=t(x)

## Creates a distance matrix. Method: Bray-Curtis, package vegan
data.dist = as.matrix(vegdist(x, "bray", upper = FALSE, na.rm = T))

## Check that the order is the same than with the metadata
table(fullMetadata_OFFS$SampleID == rownames(data.dist))

# We use a PERMANOVA to test the hypothesis that paternal treatment, 
# offspring treatment, sex and their interactions significantly influencing global methylation
perm <- how(nperm = 1000) # 1000 permutations
setBlocks(perm) <- with(fullMetadata_OFFS, Family) # define the permutation structure considering family

## Full model
adonis2(data.dist ~ PAT * outcome * Sex, data = fullMetadata_OFFS, permutations = perm)
## only PAT significant: 
#     Df  SumOfSqs  R2      F    Pr(>F) 
# PAT 1   0.04012 0.01726 1.9072 0.000999 ***

## remove the non significant interactions
adonis2(data.dist ~ PAT + outcome + Sex, data = fullMetadata_OFFS, permutations = perm)
## again only PAT significant: cluster by paternal treatment
#     Df  SumOfSqs  R2      F    Pr(>F) 
# PAT 1  0.04012 0.01726 1.9139 0.000999 ***

##### NMDS
# find the best number of dimensions (goeveg lib)
## Clarke 1993 suggests the following guidelines for acceptable stress values: <0.05 = excellent, <0.10
# = good, <0.20 = usable, >0.20 = not acceptable. The plot shows the border of the 0.20 stress value
# limit. Solutions with higher stress values should be interpreted with caution and those with stress
# above 0.30 are highly suspect
dimcheckMDS(
  data.dist,
  distance = "bray",
  k = 7,
  trymax = 100,
  autotransform = TRUE
)
abline(h = 0.1, col = "darkgreen")

#Create NMDS based on bray-curtis distances - metaMDS finds the
# most stable NMDS solution by randomly starting from different points in your data
set.seed(1234)

# generate CpG sums these will be NA is any are missing
csum <- colSums(x)
# check if any are missing
any(is.na(csum))
# TRUE
# yes, some missing, so which ones?
length(which(!is.na(csum)))#179 complete positions

x = x[,which(!is.na(csum))]

NMDS <- metaMDS(comm = x, distance = "bray", maxit=1000, k = 4)

#check to see stress of NMDS
mystressplot <- stressplot(NMDS) 

#extract plotting coordinates
MDS1 = NMDS$points[,1] ; MDS2 = NMDS$points[,2] ; MDS3 = NMDS$points[,3]
## OR #extract NMDS scores (x and y coordinates)
## data.scores = as.data.frame(scores(NMDS))

#create new data table (important for later hulls finding)
# with plotting coordinates and variables to test (dim 1,2,3)
NMDS_dt = data.table::data.table(MDS1 = MDS1, MDS2 = MDS2, MDS3 = MDS3,
                                 trtG1G2 = fullMetadata_OFFS$trtG1G2)

# generating convex hulls splitted by myvar in my metadata:
hulls <- NMDS_dt[, .SD[chull(MDS1, MDS2)], by = trtG1G2]

## to insert image (expe design)
img <- readPNG("../../data/designExpeSimple.png")
g <- rasterGrob(img, interpolate=TRUE)

myNMDSplot <- ggplot(NMDS_dt, aes(x=MDS1, y=MDS2)) +
  geom_polygon(data = hulls, aes(fill=trtG1G2), alpha=0.3) +
  # scale_color_manual(values = colOffs)+
  scale_fill_manual(values = colOffs)+
  geom_point(aes(fill=trtG1G2, shape=trtG1G2), col = "black", size = 3, alpha =.5) +
  scale_shape_manual(values = c(21,22,23,24)) +
  # geom_label(aes(label=rownames(NMDS$points), col = fullMetadata_OFFS$Family))+
  theme_bw() +
  theme(legend.position = "none") +
  annotation_custom(g, xmin=-0.25, xmax=-0.1, ymin=0.05, ymax=0.15) +
  ggtitle("NMDS analysis of the 179 CpG sites from parental control/infected DMS,\npresent in all offspring")
myNMDSplot

################################################
################ Venn diagrams #################
################################################
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
png(file="Rfigures/VennDMSinhalffish_intersectingCpGs_hyper_hypo.png", width = 5.5, height = 4.5, units = 'in', res = 300)
pushViewport(plotViewport(layout=grid.layout(2, 3)))
pushViewport(plotViewport(layout.pos.col=2, layout.pos.row=1))
grid.draw(allVenn)
popViewport()
pushViewport(plotViewport(layout.pos.col=1, layout.pos.row=2))
grid.draw(hypoVenn)
popViewport()
pushViewport(plotViewport(layout.pos.col=3, layout.pos.row=2))
grid.draw(hyperVenn)
dev.off()

rm(allVenn, hypoVenn, hyperVenn)

######################
## Features Annotation (use package genomation v1.24.0)
## NB Promoters are defined by options at genomation::readTranscriptFeatures function. 
## The default option is to take -1000,+1000bp around the TSS and you can change that. 
## -> following Heckwolf 2020 and Sagonas 2020, we consider 1500bp upstream and 500 bp downstream

#########################
## load genome annotation
gene.obj=readTranscriptFeatures("../../gitignore/bigdata/Gy_allnoM_rd3.maker_apocrita.noseq_corrected.bed12",
                                remove.unusual = FALSE, up.flank = 1500, down.flank = 500)

gene.obj.2=readTranscriptFeatures("../../gitignore/bigdata/Gy_allnoM_rd3.maker_apocrita.noseq_corrected_afterAGATcuration.bed12",
                                  remove.unusual = FALSE, up.flank = 1500, down.flank = 500)


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

DMS15pc_G1_half = myAnnotateDMS(DMS15pc_G1_half, as.data.frame(diffAnn_PAR@members))
DMS15pc_G1_half_HYPO = myAnnotateDMS(DMS15pc_G1_half[DMS_info_G1$direction %in% "hypo",],
                                     as.data.frame(myannot$G1hypo@members))
DMS15pc_G1_half_HYPER = myAnnotateDMS(DMS15pc_G1_half[DMS_info_G1$direction %in% "hyper",],
                                      as.data.frame(myannot$G1hyper@members))

DMS15pc_G2_controlG1_half = myAnnotateDMS(DMS15pc_G2_controlG1_half, as.data.frame(diffAnn_G2_controlG1@members))
DMS15pc_G2_controlG1_half_HYPO = myAnnotateDMS(DMS15pc_G2_controlG1_half[DMS_info_G2_G1c_final$direction %in% "hypo",],
                                               as.data.frame(myannot$G2G1chypo@members))
DMS15pc_G2_controlG1_half_HYPER = myAnnotateDMS(DMS15pc_G2_controlG1_half[DMS_info_G2_G1c_final$direction %in% "hyper",],
                                                as.data.frame(myannot$G2G1chyper@members))

DMS15pc_G2_infectedG1_half = myAnnotateDMS(DMS15pc_G2_infectedG1_half, as.data.frame(diffAnn_G2_infectedG1@members))
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

###############################
## Gene association with DMS ##
###############################

## Keep CpG on genes:

## Heckwolf 2020: To be associated to a gene, the pop-DMS had to be either inside the gene or, 
## if intergenic, not further than 10 kb away from the TSS. 
rows2rm <- which(diffAnn_PAR@dist.to.TSS$dist.to.feature>10000 & 
                   rowSums(diffAnn_PAR@members) %in% 0)

DMS15pc_G1_half_GENES <- DMS15pc_G1_half[-rows2rm,]

# 4393 to keep, 681 to rm
############ Follow up
gene.obj$promoters

DMS15pc_G1_half_GENES





setwd("/home/user/R/x86_64-pc-linux-gnu-library/4.0/genomation/extdata")
gff.file=system.file("extdata/Bubalus_bubalis.gtf", package = "genomation")
gff = gffToGRanges(gff.file)
head(gff)
Anno <- GRangesList(gff)
setwd("~/Desktop/Final_sort")
diffAnnhyper3=annotateWithGeneParts(as(myDiff25p.hyper,"Granges"),Anno)


?readTranscriptFeatures
head(getAssociationWithTSS(diffAnn_PAR))
