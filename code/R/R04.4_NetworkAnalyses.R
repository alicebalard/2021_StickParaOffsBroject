## Network analyses
## A. Balard
## February 2022

machine="mythinkpad" # define the machine we work on
loadALL = FALSE # only load CpG shared by half fish per trt group
loadannot = TRUE # load genome annotations
source("R02.3_DATALOAD.R")

#### NB: add BCI (was calculated only in R03.1)
## Kaufmann et al. 2014: Body condition of the G2 fish, an estimate of fish health and a predictor
# of energy reserves and reproductive success, was calculated using there residuals from the 
# regression of body mass on body length (Chellappaet al.1995).
fullMetadata_OFFS$BCI <- residuals(lmer(Wnettofin ~ Slfin * Sex + (1|Family), data=fullMetadata_OFFS))
## and for parents (no sex difference, only males):
fullMetadata_PAR$BCI <- residuals(lmer(Wnettofin ~ Slfin + (1|Family), data=fullMetadata_PAR))

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

###########################################
#### Co-methylation network construction
enableWGCNAThreads()

## Prepare datasets needed. We need the individual beta values for the DMS positions:
prepareNetworkData <- function(myMethylDiff, originalUniteCov, fullMetadata){
  ## Keep full methylation data for (1) positions only at DMS and (2) correct samples
  myMethylKit_DMS = originalUniteCov[paste(originalUniteCov$chr, originalUniteCov$start, originalUniteCov$end) %in% 
                               paste(myMethylDiff$chr, myMethylDiff$start, myMethylDiff$end)]
  myMethylKit_DMS = reorganize(methylObj = myMethylKit_DMS, 
                               treatment = myMethylDiff@treatment, 
                               sample.ids = myMethylDiff@sample.ids)
  ## Keep correct metadata
  myMetadata = fullMetadata[fullMetadata$trtG1G2_NUM %in% myMethylDiff@treatment,]
  return(list(myMethylDiff=myMethylDiff, myMethylKit_DMS = myMethylKit_DMS, myMetadata = myMetadata))
}

myG1list <- prepareNetworkData(myMethylDiff = DMS15pc_G1_ALL,
                               originalUniteCov = uniteCov6_G1_woSexAndUnknowChrOVERLAP,
                               fullMetadata = fullMetadata_PAR)

myG2_G1control_list <- prepareNetworkData(myMethylDiff = DMS15pc_G2_controlG1_ALL,
                                          originalUniteCov = uniteCov14_G2_woSexAndUnknowChrOVERLAP,
                               fullMetadata = fullMetadata_OFFS)

myG2_G1infected_list <- prepareNetworkData(myMethylDiff = DMS15pc_G2_infectedG1_ALL,
                                           originalUniteCov = uniteCov14_G2_woSexAndUnknowChrOVERLAP,
                               fullMetadata = fullMetadata_OFFS)

## Remove all but my 3 lists to avoid confusion:
rm(list = ls()[!ls() %in% c("myG1list", "myG2_G1control_list", "myG2_G1infected_list")])

#### PB: how to deal with missing data?
## https://elifesciences.org/articles/59201#s4
## "Missing values were imputed using a kNN sliding window; missing methylation values were assigned the average value of the five nearest neighbors by Euclidean distance within a 3Mb window"

## OR: network by GENE, with average methylation value of the gene.

## OR: we focus now on the CpG positions covered in ALL fish


ARG = myG2_G1infected_list
#### MAKEFUNTION HERE

### Prepare data
## https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-01-dataInput.pdf
## 1.a Loading methylation data

## Make a dataframe in which each column is a sample and each row a CpG position, 
## with values = beta value (methylation frequency = numCs/coverage)
myDF_betaVal = methylKit::getData(ARG$myMethylKit_DMS)[grep("numCs", names(ARG$myMethylKit_DMS))]/
  methylKit::getData(ARG$myMethylKit_DMS)[grep("coverage", names(ARG$myMethylKit_DMS))]
names(myDF_betaVal)= ARG$myMethylKit_DMS@sample.ids

## Transpose
myDF_betaVal=data.frame(t(myDF_betaVal))
colnames(myDF_betaVal)=paste0("CpG", colnames(myDF_betaVal))

## 1.b Checking data for excessive missing values and outlier samples
gsg = goodSamplesGenes(myDF_betaVal, verbose = 3, minFraction = 0.9);
gsg$allOK

# If the last statement returns TRUE, all genes have passed the cuts. If not, 
# we remove the offending genes and samples from the data
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(myDF_betaVal)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(myDF_betaVal)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  myDF_betaVal = myDF_betaVal[gsg$goodSamples, gsg$goodGenes]
}

## Cluster the samples (in contrast to clustering genes that will come later) to see if there are any obvious
## outliers:
sampleTree = hclust(dist(myDF_betaVal), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2) ## All good!

## 1.c Loading clinical trait data
## Form a data frame analogous to methylation data that will hold the individuals traits
## Resistance: load of parasites: myMetaData$No.Worms
## Tolerance: BCI 
ARG$myMetadata$Family_NUM = as.numeric(gsub("Fam", "",ARG$myMetadata$Family))
traitRows = match(rownames(myDF_betaVal), ARG$myMetadata$SampleID);
datTraits = (ARG$myMetadata)[traitRows,];
rownames(datTraits) = ARG$myMetadata[traitRows, "ID"];
datTraits = datTraits[c("trtG1G2_NUM", "Family_NUM", "No.Worms", "BCI")] ## here keep interesting CONTINUOUS traits
collectGarbage();

## We now have the methylation data in the variable datMeth, and the corresponding traits in the variable
## datTraits. Before we continue with network construction and module detection, we visualize how
## the traits relate to the sample dendrogram
# Re-cluster samples
sampleTree2 = hclust(dist(myDF_betaVal), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap\nClustering dendrogram of samples based on their Euclidean distance.")

## 2.Automatic network construction and module detection
## https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-02-networkConstr-auto.pdf

## 2.1 Choosing the soft-thresholding power: analysis of network topology
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(myDF_betaVal, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5) ; par(mfrow = c(1,2)); cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.9,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
abline(h=30,col="red")
abline(h=100,col="red")
## Figure 1: Analysis of network topology for various soft-thresholding powers. The left panel shows the scale-free fit
# index (y-axis) as a function of the soft-thresholding power (x-axis). The right pabnel displays the mean connectivity
# (degree, y-axis) as a function of the soft-thresholding power (x-axis)

## 2.2 One-step network construction and module detection
net = blockwiseModules(myDF_betaVal, power = 3, # 12 is default for signed network
                       replaceMissingAdjacencies = TRUE,
                       TOMType = "signed", # preserved direction +/- of correlations
                       saveTOMs = TRUE, checkMissingData = TRUE, 
                       corType = "bicor", # handles outliers better than Pearson
                       saveTOMFileBase = "G1TOM",
                       verbose = 3)

# We now return to the network analysis. To see how many modules were identified and what the module sizes are,
# one can use table(net$colors). Its output is
table(net$colors) # The label grey is reserved for CpGs outside of all modules

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
# 
################ GRAPHS ###################
library(tidyverse); library(igraph); library(ggraph); library(tidygraph); library(ggregplot)
library(fs); library(magrittr); library(cowplot); library(ggplot2)
library(patchwork); library(Matrix); library(MASS); library(reshape2)

## Export from WGCNA to igraph for visualisation

# Recalculate topological overlap
TOM = TOMsimilarityFromExpr(myDF_betaVal, power = 12);

rownames(TOM) <- names(myDF_betaVal)
colnames(TOM) <- names(myDF_betaVal)
# this is an adjacency matrix

## Make graph:
myGraph <- TOM %>%
  # stores the correlations as edge weight:
  graph_from_adjacency_matrix(mode = "undirected", weighted = T) %>% # igraph. Weight saved in edges
  # as_tbl_graph(nodes = CpGInfo) %>% # %>% # tidygraph with nodes CpG info
  simplify()%>% # removes self-loops
  as_tbl_graph() %>% # tidygraph with nodes CpG name
  activate(nodes) %>%
#  left_join(CpGInfo, by = c("name" = "CpGpos"))%>% # add CpGInfo: chromosome location
  left_join(data.frame(name= names(moduleLabels), color = moduleColors), by = c("name" = "name")) # add previously detected module colors

## add layout based on algorythms, compare different layouts
Layout1 <- layout_with_fr(myGraph, # Fruchterman-Reingold as an example
                          weights = myGraph %>%
                            activate(edges) %>%
                            pull(weight)) # Adding in weights explicitly
# 
# Layout2 <- layout_with_kk(myGraph,
#                           weights = myGraph %>%
#                             activate(edges) %>%
#                             pull(weight)) # Kamada-Kawai
# Layout3 <- layout_with_drl(myGraph,
#                            weights = myGraph %>%
#                              activate(edges) %>%
#                              pull(weight)) # Deep reinforcement learning
# Layout4 <- layout_with_mds(myGraph) # Multi-dimensional scaling
# 
# # Comparing the layout algorithms ####
# list(Layout1, Layout2, Layout3, Layout4) %>%
#   map(~myGraph %>% # Going through the list of layouts
#         ggraph(.x) + # Applying a new layout to the graph each time
#         geom_edge_link0(alpha = 0.1) +
#         geom_node_point() +
#         theme_bw() +
#         coord_fixed()) %>%
#   ArrangeCowplot() +  # My function to plot them all together using patchwork
#   plot_layout(ncol = 2) + # 2 columns
#   plot_annotation(tag_levels = "A") # Label them with letters

## Plot:
myGraph %>% # my igraph - tidygraph
  ggraph(Layout1) +  # then ggraph
  geom_edge_link(aes(alpha=weight)) +
  scale_edge_alpha(range=c(0.1,1))+
  geom_node_point(aes(color=color), size = 7) +
  #  scale_colour_manual(values = c("grey", "turquoise"))+
  geom_node_text(aes(label = name), size = 3)


# ########### Network "bonus" ideas
# 
# ## test other nodes characteristics
# NodeTraitGet <- function(Graph, WeightVar = "weight"){
#   Weights <- Graph %>% activate(edges) %>% as.data.frame
#   Weights <- Weights[,WeightVar]
#   Graph %>% activate(edges) %>% mutate(weight = Weights)
#   AM <- Graph %>% 
#     get.adjacency(attr = WeightVar) %>% 
#     as.matrix
#   list(
#     ID = colnames(AM),
#     Degree = colSums(AM > 0),
#     Strength = colSums(AM),
#     Eigenvector = Graph %>% eigen_centrality(weights = NULL) %>% extract2("vector"),
#     Eigenvector_Weighted = Graph %>% 
#       eigen_centrality(weights = Graph %>% activate(edges) %>% pull(weight)) %>% 
#       extract2("vector"),
#     Clustering = Graph %>% transitivity(type = "local"),
#     Betweenness = Graph %>% betweenness,
#     Closeness = Graph %>% closeness
#   )
# }
# 
# Nodes <- myGraph %>% activate(nodes) %>% pull(name)
# 
# myGraph %<>% NodeTraitGet() %>% bind_cols %>% 
#   mutate(name = Nodes) %>% 
#   left_join(
#     myGraph %>% activate(nodes), .,
#     by = "name"
#   )
# 
# ## Color according to e.g. degree of connection to the other nodes
# ## Plot:
# myGraph %>% # my igraph - tidygraph
#   ggraph(Layout1) +  # then ggraph
#   geom_edge_link(aes(alpha=weight)) + 
#   scale_edge_alpha(range=c(0.1,1))+
#   geom_node_point(aes(color=Degree), size = 5) +
#   scale_colour_gradient(low = "blue", high = "red")+
#   geom_node_text(aes(label = name), size = 3) 
# 
# ################ WGCNA
