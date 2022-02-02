## Differential methylation analyses
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
## Metylation metadata merged with Kostas files
fullMetadata <- read.csv("../../data/fullMetadata137_Alice.csv")
## relevel treatments for graphs
fullMetadata$trtG1G2 <- factor(as.factor(fullMetadata$trtG1G2), levels = c("Control", "Exposed","NE_control", "NE_exposed", "E_control", "E_exposed"  ))
fullMetadata$trtG1G2_NUM <- as.numeric(as.factor(fullMetadata$trtG1G2))
## family as factor for models
fullMetadata$Family <- as.factor(fullMetadata$Family)
# paternal exposure
fullMetadata$PAT="Exposed father group"
fullMetadata$PAT[fullMetadata$trtG1G2 %in% c("Control", "NE_control", "NE_exposed")]="Control father group"

## only parents
fullMetadata_PAR <- fullMetadata[fullMetadata$Generat %in% "P",]
## only offspring
fullMetadata_OFFS <- fullMetadata[fullMetadata$Generat %in% "O",]
# without sample 22 outlier that block the view on PCA
fullMetadata_OFFS_no22 <- fullMetadata_OFFS[!fullMetadata_OFFS$ID %in% "S22",]

## Load previously united data (all 6 treatments)
## uniteCovALL: CpG covered in ALL individuals (has no NAs, useful for exploratory clustering analyses)
#load("/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/Data/05MethylKit/output/uniteCovALL_woSexAndUnknownChr.RData")
base::load("../../gitignore/output/uniteCovALL_woSexAndUnknownChr.RData")
## rename
fullMethylKitObj = uniteCovALL_woSexAndUnknowChr; rm(uniteCovALL_woSexAndUnknowChr)

## For further analyses: CpG covered in at least 2 individuals per group 
## (Kostas took 2; Melanie is more stringent; let's see which makes sense)
load("../../gitignore/output/uniteCov2_woSexAndUnknownChr.RData")
## rename
AL2MethylKitObj = uniteCov2_woSexAndUnknowChr ; rm(uniteCov2_woSexAndUnknowChr)

#load("/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/Data/05MethylKit/output/uniteCov2_woSexAndUnknownChr.RData")
#load("/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/Data/05MethylKit/output/uniteCov6_woSexAndUnknowChr.RData")

##### Create methylkit objects for different analyses:
# create a methylKit object with ONLY the parents - positions shared by AT LEAST 2 INDIVIDUALS
uniteCov2_woSexAndUnknowChr_PAR=reorganize(
  AL2MethylKitObj,
  sample.ids=fullMetadata_PAR$ID,
  treatment=fullMetadata_PAR$trtG1G2_NUM)

## Remove positions with only NA
# methylBase object with 1772152 rows
test <- methylKit::getData(uniteCov2_woSexAndUnknowChr_PAR)
A <- rowSums(is.na(test[grep("coverage", names(test))]))
table(A)#600226 columns have only NAs! Remove that
mynonNAcols <- which(rowSums(is.na(test[grep("coverage", names(test))]))!=24)
## SUBSET TO REMOVE FULL NA ROWS!!!!:
uniteCov2_woSexAndUnknowChr_PAR <- uniteCov2_woSexAndUnknowChr_PAR[mynonNAcols]
rm(test, A,mynonNAcols)

# create a methylKit object with ONLY the offspring - positions shared by ALL
uniteCovALL_woSexAndUnknowChr_OFF=reorganize(
  fullMethylKitObj,
  sample.ids=fullMetadata_OFFS$ID,
  treatment=fullMetadata_OFFS$trtG1G2_NUM)

# create a methylKit object with ONLY the offspring - positions shared by AT LEAST 2 INDIVIDUALS
uniteCov2_woSexAndUnknowChr_OFF=reorganize(
  AL2MethylKitObj,
  sample.ids=fullMetadata_OFFS$ID,
  treatment=fullMetadata_OFFS$trtG1G2_NUM)

# without sample 22 outlier that block the view on PCA
uniteCovALL_woSexAndUnknowChr_OFF_no22=reorganize(
  uniteCovALL_woSexAndUnknowChr_OFF,
  sample.ids=fullMetadata_OFFS_no22$ID,
  treatment=fullMetadata_OFFS_no22$trtG1G2)

##############################
##### Analysis workflow: #####
# PART 1.Methylation profile: mixed model to test paternal + family effect;
#   or Adonis multivariate stats, presence/abs, block effect=family
# PART 2.DMS
# PART 3.DMR: Precise localisation of blocks investing in methylation
# PART 4.Network: Explore other parameters than geographic: are the modules
#   on a high recombination place? Or with high mutation rate? 
#   What does modules capture on top of DMRs?
##### Then 
# - Gene Ontology
# - Final: Check robustness: reshuffle the controls and check that we find less DMS/DMR/modules
# 
# These analyses will be done on:
# A. Parents ctr-trt
# B. Offspring ctr-trt pat1 & ctrl-trt pat2 (different offspring treatments)
# C. Offspring ctrl-ctrl & trt-trt (different paternal treatment)
##############################

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
# pdf("../../data/fig/clusterALLCpG_offspings.pdf", width = 17, height = 7)
makePrettyMethCluster(uniteCovALL_woSexAndUnknowChr_OFF, fullMetadata_OFFS,
                      my.cols.trt=c("#ffe680ff","#ff6600ff", "#aaccffff", "#aa00d4ff"),
                      my.cols.fam = c(1:4))
# dev.off()

############################
## HERE CHOOSE which dataset
dataset = uniteCovALL_woSexAndUnknowChr_OFF; metadata = fullMetadata_OFFS

# creates a matrix containing percent methylation values NB rm sample 22!!
perc.meth=percMethylation(dataset)
# KOSTAS MBE: "Methylated sites and regions with low variation
# and a standard deviation below 0.3, that is, noninformative
# sites across individuals, were excluded from the cluster analyses"
SD=apply(perc.meth,1, sd, na.rm = TRUE)
perc.meth <- perc.meth[-which(SD<0.3),]

x=t(perc.meth)

# creates a distance matrix. Method: Bray-Curtis, package vegan
data.dist = as.matrix((vegdist(x, "bray", upper = FALSE))) 

# We use a PERMANOVA to test the hypothesis that paternal treatment,
# family and sex induced changes in genome-wide methylation. 
# We used Euclidian distance matrix between methylation frequencies
# and the function Adonis from the R package vegan:

## Within each family, are paternal treatment, offspring treatment, sex and their interactions
## significantly influencing global methylation?
perm <- how(nperm = 1000) # 1000 permutations
setBlocks(perm) <- with(metadata, Family) # define the permutation structure

## Full model
adonis2(data.dist ~ PAT * outcome * Sex, data = metadata, permutations = perm)

## remove the non significant interactions
adonis2(data.dist ~ PAT + outcome + Sex, data = metadata, permutations = perm)

## Remove 1 factor by turn - backwards simplification
adonis2(data.dist ~ outcome + Sex, data = metadata, permutations = perm)

adonis2(data.dist ~ PAT + Sex, data = metadata, permutations = perm)

adonis2(data.dist ~ PAT + outcome, data = metadata, permutations = perm)

## --> We found significant differences in global methylation due to
## paternal treatment and sex; outcome is not. 

########## NMDS
# find the best number of dimensions
library(goeveg)
## Clarke 1993 suggests the following guidelines for acceptable stress values: <0.05 = excellent, <0.10
# = good, <0.20 = usable, >0.20 = not acceptable. The plot shows the border of the 0.20 stress value
# limit. Solutions with higher stress values should be interpreted with caution and those with stress
# above 0.30 are highly suspect
# dimcheckMDS(
#   data.dist,
#   distance = "bray",
#   k = 7,
#   trymax = 100,
#   autotransform = TRUE
# )
# abline(h = 0.1, col = "darkgreen")
# Goodness of fit for NMDS suggested the presence of six dimensions
# with a stress value <0.1

#Create NMDS based on bray-curtis distances - metaMDS finds the
# most stable NMDS solution by randomly starting from different points in your data
set.seed(123)
NMDS <- metaMDS(comm = x, distance = "bray",maxit=1000, k = 6)

#check to see stress of NMDS
stressplot(NMDS) # very not linear when ID 22 added!!

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

## with paternal treatment
find_hull <- function(my_data) my_data[chull(my_data[,1], my_data[,2]), ]
#function to create convex hulls, though ggplot has the stat_ellipse() function that can do this automatically.
hulls <- ddply(NMDS2, "PAT", find_hull)
#generating convex hulls by "PAT" in my metadata
ggplot(NMDS2, aes(x=MDS1, y=MDS2)) +
  geom_polygon(data = hulls, aes(col=PAT), alpha=0) +
  scale_color_manual(values = c("grey","yellow"))+
  scale_fill_manual(values = c("grey","yellow"))+
  geom_point(aes(fill=PAT, shape=outcome), size = 3, alpha = .6) +
  geom_label(aes(label=row.names(NMDS2)))+
  scale_shape_manual(values = c(21,22)) +
  theme_bw() +
  theme(legend.title=element_blank(), legend.position = "top")

## with Sex
hulls2 <- ddply(NMDS2, "Sex", find_hull)
ggplot(NMDS2, aes(x=MDS1, y=MDS2)) +
  geom_polygon(data = hulls2, aes(col=Sex), alpha=0) +
  scale_color_manual(values = c("red","blue"))+
  scale_fill_manual(values = c("red","blue"))+
  geom_point(aes(fill=Sex, shape=outcome), size = 3, alpha = .6) +
  scale_shape_manual(values = c(21,22)) +
  theme_bw() +
  theme(legend.title=element_blank(), legend.position = "top")

###############################################
### PART 2: Differential methylation sites ####
###############################################
options(stringsAsFactors = FALSE);
table(fullMetadata_OFFS$trtG1G2,fullMetadata_OFFS$trtG1G2_NUM)

##########################################################
## Comparison 1: BASELINE -> Parents (control vs infected) 
## Extract DMS from parents (at least in 2 fish), annotate them, compare with Kostas results

## Calculate DMS accounting for covariates: family
cov = data.frame(Family = fullMetadata_PAR$Family)
myDiffMeth=calculateDiffMeth(uniteCov2_woSexAndUnknowChr_PAR, covariates = cov, mc.cores = 4)

# We select the bases that have q-value<0.01 and percent methylation difference larger than 15%.
# NB: arg type="hyper" or type="hypo" gives hyper-methylated or hypo-methylated regions/bases.
myDiff1_15p = getMethylDiff(myDiffMeth,difference=15,qvalue=0.01)

myDiff1_15p # 6544 positions
# saveRDS(myDiff1_15p, file = "../../gitignore/output/myDiff1_15p_parentalDiffMeth.RDS")

myDiff1_15p <- readRDS("../../gitignore/output/myDiff1_15p_parentalDiffMeth.RDS")

# annotation!!
library(genomation)

#library(GenomicFeatures)

gene.obj=readTranscriptFeatures("../../gitignore/bigdata/Gy_allnoM_rd3.maker_apocrita.noseq_corrected.bed12", remove.unusual = FALSE)

# annotate differentially methylated CpGs with promoter/exon/intron using annotation data
annotateWithGeneParts(as(myDiff1_15p,"GRanges"),gene.obj)

## Kostas MBE: The DMSs and regions were predominately
#found in intergenic regions (47.74% and 48.94%, respecti vely),with introns (26.19% and 23.09), exons (15.07% and 13.98%),and promoters (11% and 13.98%) showing lower proportions

?annotateWithGeneParts()


### YOU'RE HERE ;)

###################################################################################################
## Comparison 2: offsprings infected from unifected (trtgroup 4) and infected (trt group 6) fathers
## NB CpG present in ALL FISH (not at least 2) for a start

### Select metadata
myMetaData = fullMetadata_OFFS[fullMetadata_OFFS$trtG1G2_NUM %in% c(4,6),]

### Select methylKit object
myMethylKit = reorganize(fullMethylKitObj, 
                         sample.ids=myMetaData$ID,
                         treatment=myMetaData$trtG1G2_NUM)

### Calculate DMS accounting for covariates: family
cov = data.frame(Family = myMetaData$Family)
myDiffMeth=calculateDiffMeth(myMethylKit, covariates = cov, mc.cores = 4)

# We select the bases that have q-value<0.01 and percent methylation difference larger than 15%.
# NB: arg type="hyper" or type="hypo" gives hyper-methylated or hypo-methylated regions/bases.
myDiff2_15p = getMethylDiff(myDiffMeth,difference=15,qvalue=0.01)

myDiff2_15p # 66 positions

##########################################################
#### Annotating differentially methylated bases or regions
#BiocManager::install("genomation")
# library(genomation)
## TBC

##########################################################
#### Co-methylation network construction
### Help: https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-01-dataInput.pdf
library(tidyverse)     # tidyverse will pull in ggplot2, readr, other useful libraries
library(magrittr)      # provides the %>% operator
#BiocManager::install("WGCNA") 
library(WGCNA) 

#######
## 1.Data preparation

## Extract the targeted CpG positions from the methylkit object
myMethylKit$position <- paste(myMethylKit$chr, myMethylKit$start, myMethylKit$end)
myDF <- myMethylKit[myMethylKit$position %in% paste(myDiff1_15p$chr, myDiff1_15p$start, myDiff1_15p$end),]

## Make myDF_fr : a dataframe in which each column is a sample or auxilliary info,
## and each row a CpG position, with values = beta value (methylation frequency = numCs/coverage)
myDF <-methylKit::getData(myDF)
myDF_fr <- data.frame(matrix(ncol=length(myMethylKit@sample.ids), nrow=nrow(myDF)))
names(myDF_fr) = myMethylKit@sample.ids
for(i in 1:ncol(myDF_fr)){
  myDF_fr[i] <- myDF[paste0("numCs", i)]/myDF[paste0("coverage", i)]
}

## Store CpG information (Chr, start, end)
CpGInfo <- data.frame(CpGpos = paste0("CpG", 1:length(myDF$chr)), 
                      Chromosome = myDF$chr,
                      start = myDF$start,
                      end=myDF$end)

## We now transpose the expression data for further analysis.
datMeth = as.data.frame(t(myDF_fr));
names(datMeth) = CpGInfo$CpGpos

## We first check for genes and samples with too many missing values:
gsg = goodSamplesGenes(datMeth, verbose = 3);
gsg$allOK # all good!

## Cluster the samples (in contrast to clustering genes that will come later) to see if there are any obvious
## outliers:
sampleTree = hclust(dist(datMeth), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2) ## All good!

## Form a data frame analogous to expression data that will hold the individuals traits

nrow(myMetaData)

traitRows = match(rownames(datMeth), myMetaData$SampleID);
datTraits = myMetaData[traitRows,];
rownames(datTraits) = myMetaData[traitRows, "ID"];
datTraits = datTraits[c("BCI", "integration.RLU.")] ## here keep interesting CONTINUOUS traits (no "Family", "Sex", "trtG1G2")
collectGarbage();

## We now have the expression data in the variable datMeth, and the corresponding traits in the variable
## datTraits. Before we continue with network construction and module detection, we visualize how 
## the traits relate to the sample dendrogram
# Re-cluster samples
sampleTree2 = hclust(dist(datMeth), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")

## if needed
## save(datMeth, datTraits, file = "filename.RData")

#######
## 2.Automatic network construction and module detection

###
## 2.1 Choosing the soft-thresholding power: analysis of network topology
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datMeth, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5) ; par(mfrow = c(1,2)); cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
## Figure 1: Analysis of network topology for various soft-thresholding powers. The left panel shows the scale-free fit
# index (y-axis) as a function of the soft-thresholding power (x-axis). The right panel displays the mean connectivity
# (degree, y-axis) as a function of the soft-thresholding power (x-axis)
## We choose the power 3, which is the lowest power for which the scale-free topology fit
# index curve flattens out upon reaching a high value.

###
## 2.2 One-step network construction and module detection
net = blockwiseModules(datMeth, power = 3,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "femaleMouseTOM",
                       verbose = 3)

# We now return to the network analysis. To see how many modules were identified and what the module sizes are,
# one can use table(net$colors). Its output is
table(net$colors) # 1 module, with 42 CpG. The label 0 is reserved for CpGs outside of all modules

# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
# save(MEs, moduleLabels, moduleColors, geneTree,
#      file = "FemaleLiver-02-networkConstruction-auto.RData")

################ GRAPHS ###################
library(tidyverse); library(igraph); library(ggraph); library(tidygraph); library(ggregplot)
library(fs); library(magrittr); library(cowplot); library(ggplot2)
library(patchwork); library(Matrix); library(MASS); library(reshape2)

## Export from WGCNA to igraph for visualisation

# Recalculate topological overlap
TOM = TOMsimilarityFromExpr(datMeth, power = 3);

rownames(TOM) <- names(datMeth)
colnames(TOM) <- names(datMeth)
## 66 by 66: this is an adjacency matrix

## Make graph:
myGraph <- TOM %>% 
  # stores the correlations as edge weight:
  graph_from_adjacency_matrix(mode = "undirected", weighted = T) %>% # igraph. Weight saved in edges
  # as_tbl_graph(nodes = CpGInfo) %>% # %>% # tidygraph with nodes CpG info
  simplify()%>% # removes self-loops
  as_tbl_graph() %>% # tidygraph with nodes CpG name
  activate(nodes) %>% 
  left_join(CpGInfo, by = c("name" = "CpGpos"))%>% # add CpGInfo: chromosome location 
  left_join(data.frame(name= names(moduleLabels), color = moduleColors), by = c("name" = "name")) # add previously detected module colors

## add layout based on algorythms, compare different layouts
Layout1 <- layout_with_fr(myGraph, # Fruchterman-Reingold as an example
                          weights = myGraph %>% 
                            activate(edges) %>% 
                            pull(weight)) # Adding in weights explicitly

Layout2 <- layout_with_kk(myGraph, 
                          weights = myGraph %>% 
                            activate(edges) %>% 
                            pull(weight)) # Kamada-Kawai
Layout3 <- layout_with_drl(myGraph, 
                           weights = myGraph %>% 
                             activate(edges) %>% 
                             pull(weight)) # Deep reinforcement learning
Layout4 <- layout_with_mds(myGraph) # Multi-dimensional scaling

# Comparing the layout algorithms ####
list(Layout1, Layout2, Layout3, Layout4) %>% 
  map(~myGraph %>% # Going through the list of layouts
        ggraph(.x) + # Applying a new layout to the graph each time
        geom_edge_link0(alpha = 0.1) +
        geom_node_point() +
        theme_bw() +
        coord_fixed()) %>% 
  ArrangeCowplot() +  # My function to plot them all together using patchwork
  plot_layout(ncol = 2) + # 2 columns
  plot_annotation(tag_levels = "A") # Label them with letters

## Plot:
myGraph %>% # my igraph - tidygraph
  ggraph(Layout1) +  # then ggraph
  geom_edge_link(aes(alpha=weight)) + 
  scale_edge_alpha(range=c(0.1,1))+
  geom_node_point(aes(color=color), size = 7) +
  scale_colour_manual(values = c("grey", "turquoise"))+
  geom_node_text(aes(label = name), size = 3) 


########### Network "bonus" ideas

## test other nodes characteristics
NodeTraitGet <- function(Graph, WeightVar = "weight"){
  Weights <- Graph %>% activate(edges) %>% as.data.frame
  Weights <- Weights[,WeightVar]
  Graph %>% activate(edges) %>% mutate(weight = Weights)
  AM <- Graph %>% 
    get.adjacency(attr = WeightVar) %>% 
    as.matrix
  list(
    ID = colnames(AM),
    Degree = colSums(AM > 0),
    Strength = colSums(AM),
    Eigenvector = Graph %>% eigen_centrality(weights = NULL) %>% extract2("vector"),
    Eigenvector_Weighted = Graph %>% 
      eigen_centrality(weights = Graph %>% activate(edges) %>% pull(weight)) %>% 
      extract2("vector"),
    Clustering = Graph %>% transitivity(type = "local"),
    Betweenness = Graph %>% betweenness,
    Closeness = Graph %>% closeness
  )
}

Nodes <- myGraph %>% activate(nodes) %>% pull(name)

myGraph %<>% NodeTraitGet() %>% bind_cols %>% 
  mutate(name = Nodes) %>% 
  left_join(
    myGraph %>% activate(nodes), .,
    by = "name"
  )

## Color according to e.g. degree of connection to the other nodes
## Plot:
myGraph %>% # my igraph - tidygraph
  ggraph(Layout1) +  # then ggraph
  geom_edge_link(aes(alpha=weight)) + 
  scale_edge_alpha(range=c(0.1,1))+
  geom_node_point(aes(color=Degree), size = 5) +
  scale_colour_gradient(low = "blue", high = "red")+
  geom_node_text(aes(label = name), size = 3) 

################ WGCNA



