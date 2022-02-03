## Global methylation analyses
## A. Balard
## November 2021

#################### Data load & preparation ####################
library(plyr) # for join (keep row order)
library(dendextend) # The package dendextend contains many functions for changing the appearance of a dendrogram and for comparing dendrograms.
library(ggplot2)
library(ggsignif) ## for significance bars on ggplot
library(lme4) ## for mixed models
library(nlme) ## for mixed models
library(tidyverse)
library(emmeans) ## for post-hoc Tukey tests
library(methylKit)
library(vegan) ## for Adonis

## load custom functions
source("customRfunctions.R")

## Load samples metadata
source("R01.3_prepMetadata.R")

## Load methylation data
source("R01.4_prepMethyldata.R")

#############################################################
### PART 1: Methylation profiles, CpG present in all fish ###
#############################################################

## Dendogram of methylations
## All samples:
# pdf("../../data/fig/clusterALLCpG.pdf", width = 16, height = 7)
# makePrettyMethCluster(fullMethylKitObj, fullMetadata,
#                       my.cols.trt=c("#333333ff","#ff0000ff","#ffe680ff","#ff6600ff","#aaccffff","#aa00d4ff"),
#                       my.cols.fam = c(1:4))
# dev.off()

## offspring:
pdf("../../data/fig/clusterALLCpG_offspings.pdf", width = 17, height = 7)
makePrettyMethCluster(uniteCovALL_woSexAndUnknowChr_OFF, fullMetadata_OFFS,
                      my.cols.trt=c("#ffe680ff","#ff6600ff", "#aaccffff", "#aa00d4ff"),
                      my.cols.fam = c(1:4))
dev.off()

############################
makePercentMetMat <- function(dataset){
  # creates a matrix containing percent methylation values
  perc.meth=percMethylation(dataset)
  # KOSTAS MBE: "Methylated sites and regions with low variation
  # and a standard deviation below 0.3, that is, noninformative
  # sites across individuals, were excluded from the cluster analyses"
  SD=apply(perc.meth,1, sd, na.rm = TRUE)
  perc.meth <- perc.meth[-which(SD<0.3),]
  x=t(perc.meth)
  return(x)
}

makeDatadistFUN <- function(dataset){
  x=makePercentMetMat(dataset)
  # creates a distance matrix. Method: Bray-Curtis, package vegan
  data.dist = as.matrix((vegdist(x, "bray", upper = FALSE))) 
}

myadonisFUN <- function(dataset, metadata){
  # make distance matrix with B-C distances
  data.dist = makeDatadistFUN(dataset)
  
  # We use a PERMANOVA to test the hypothesis that paternal treatment, 
  # offspring treatment, sex and their interactions significantly influencing global methylation
  perm <- how(nperm = 1000) # 1000 permutations
  setBlocks(perm) <- with(metadata, Family) # define the permutation structure considering family
  
  ## Full model
  print(adonis2(data.dist ~ PAT * outcome * Sex, data = metadata, permutations = perm))
  ## remove the non significant interactions
  print(adonis2(data.dist ~ PAT + outcome + Sex, data = metadata, permutations = perm))
  ## Remove 1 factor by turn - backwards simplification
  # adonis2(data.dist ~ outcome + Sex, data = metadata, permutations = perm)
  # adonis2(data.dist ~ PAT + Sex, data = metadata, permutations = perm)
  # adonis2(data.dist ~ PAT + outcome, data = metadata, permutations = perm)
}

###################
## Let's run Adonis
myadonisFUN(dataset = uniteCovALL_woSexAndUnknowChr_OFF, metadata = fullMetadata_OFFS)

########## NMDS
myGOF.NMDS.FUN <- function(dataset){
  # make distance matrix with B-C distances
  data.dist = makeDatadistFUN(dataset)
  # find the best number of dimensions
  library(goeveg)
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
}

#### RUN Goodness of fit
# myGOF.NMDS.FUN(dataset = uniteCovALL_woSexAndUnknowChr_OFF) # Goodness of fit for NMDS 
# suggested the presence of six dimensions with a stress value <0.1

myNMDSplots <- function(dataset, metadata){
  ## make percent methylation matrix
  x=makePercentMetMat(dataset)
  
  #Create NMDS based on bray-curtis distances - metaMDS finds the
  # most stable NMDS solution by randomly starting from different points in your data
  set.seed(2)
  NMDS <- metaMDS(comm = x, distance = "bray", maxit=1000, k = 6)
  
  #check to see stress of NMDS
  mystressplot <- stressplot(NMDS) 
  
  #extract plotting coordinates
  MDS1 = NMDS$points[,1]
  MDS2 = NMDS$points[,2]
  ## OR #extract NMDS scores (x and y coordinates)
  ## data.scores = as.data.frame(scores(NMDS))
  
  #create new dataframe with plotting coordinates and variables to test
  NMDS2 = data.frame(MDS1 = MDS1, MDS2 = MDS2, ID = metadata$ID,
                     PAT=as.factor(metadata$PAT), 
                     outcome=as.factor(metadata$outcome), 
                     Sex = as.factor(metadata$Sex))
  
  #function to create convex hulls, though ggplot has the stat_ellipse() function that can do this automatically.
  find_hull <- function(my_data) my_data[chull(my_data[,1], my_data[,2]), ]
  
  ## with paternal treatment
  hulls <- ddply(NMDS2, "PAT", find_hull)
  #generating convex hulls by "PAT" in my metadata
  myPATplot <- ggplot(NMDS2, aes(x=MDS1, y=MDS2)) +
    geom_polygon(data = hulls, aes(col=PAT), alpha=0) +
    scale_color_manual(values = c("grey","yellow"))+
    scale_fill_manual(values = c("grey","yellow"))+
    geom_point(aes(fill=PAT, shape=PAT), size = 3, alpha = .6) +
    #  geom_label(aes(label=row.names(NMDS2)))+
    scale_shape_manual(values = c(21,22)) +
    theme_bw() +
    theme(legend.title=element_blank(), legend.position = "top")
  
  ## with Sex
  hulls2 <- ddply(NMDS2, "Sex", find_hull)
  mySEXplot <- ggplot(NMDS2, aes(x=MDS1, y=MDS2)) +
    geom_polygon(data = hulls2, aes(col=Sex), alpha=0) +
    scale_color_manual(values = c("red","blue"))+
    scale_fill_manual(values = c("red","blue"))+
    geom_point(aes(fill=Sex, shape=Sex), size = 3, alpha = .6) +
    scale_shape_manual(values = c(21,22)) +
    theme_bw() +
    theme(legend.title=element_blank(), legend.position = "top")
  
  ## with offspring trt
  hulls2 <- ddply(NMDS2, "outcome", find_hull)
  myOFFplot <- ggplot(NMDS2, aes(x=MDS1, y=MDS2)) +
    geom_polygon(data = hulls2, aes(col=outcome), alpha=0) +
    scale_color_manual(values = c("grey","red"))+
    scale_fill_manual(values = c("grey","red"))+
    geom_point(aes(fill=outcome, shape=outcome), size = 3, alpha = .6) +
    scale_shape_manual(values = c(21,22)) +
    theme_bw() +
    theme(legend.title=element_blank(), legend.position = "top")
  return(list(mystressplot, myPATplot, mySEXplot, myOFFplot))
}

myNMDSplots(dataset = uniteCovALL_woSexAndUnknowChr_OFF, metadata = fullMetadata_OFFS)