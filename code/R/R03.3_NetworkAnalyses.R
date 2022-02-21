## Network analyses
## A. Balard
## February 2022

#################### Data load & preparation ####################
## define in which machine we're working (apocrita or mythinkpad)
##machine="apocrita"
machine="mythinkpad"
## Load DMS script: needed for here!
source("R03.2_differentialMethylation.R")

## Prepare datasets needed, remove the rest:
prepareNetworkData <- function(DMS, uniteCov, fullMetadata){
  myMethylDiff = DMS
  ## Keep full methylation data for (1) positions only at DMS and (2) correct samples
  samples = myMethylDiff@sample.ids
  pos = myMethylDiff$pos
  trt = myMethylDiff@treatment
  myMethylKit_DMS = uniteCov[paste(uniteCov$chr, uniteCov$start, uniteCov$end) %in% pos]
  myMethylKit_DMS = reorganize(methylObj = myMethylKit_DMS, treatment = trt, sample.ids = samples)
  ## Keep correct metadata
  myMetadata = fullMetadata[fullMetadata$trtG1G2_NUM %in% trt,]
  return(list(myMethylDiff=myMethylDiff, myMethylKit_DMS = myMethylKit_DMS, myMetadata = myMetadata))
}

myG1list <- prepareNetworkData(DMS = DMS15pc_PAR_half_intersect,
                               uniteCov = uniteCov6_G1_woSexAndUnknowChr,
                               fullMetadata = fullMetadata_PAR)

myG2_G1control_list <- prepareNetworkData(DMS = DMS15pc_G2_controlG1_half_intersect,
                               uniteCov = uniteCov14_G2_woSexAndUnknowChr,
                               fullMetadata = fullMetadata_OFFS)

myG2_G1infected_list <- prepareNetworkData(DMS = DMS15pc_G2_infectedG1_half_intersect,
                               uniteCov = uniteCov14_G2_woSexAndUnknowChr,
                               fullMetadata = fullMetadata_OFFS)

## Remove all but my 3 lists to avoid confusion:
rm(list = ls()[!ls() %in% c("myG1list", "myG2_G1control_list", "myG2_G1infected_list")])



## YOU're here


##########################################################
#### Co-methylation network construction
### Help: https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-01-dataInput.pdf
enableWGCNAThreads()






################# Previous

# 
# myMethylKit = uniteCov2_woSexAndUnknowChr_PAR
# # load file with your DMS
# DMS_PAR <- readRDS("../../data/myDiff1_15p_parentalDiffMeth.RDS")
# myDMS = DMS_PAR # 13737 observation
# myMetaData <- fullMetadata_PAR

## in step gsg$allOK: Removing samples: S70, S95


#######
## 1.Data preparation

######################################
## METHYLATION DATA

## Extract the targeted CpG positions from the methylkit object
myMethylKit$position <- paste(myMethylKit$chr, myMethylKit$start, myMethylKit$end)
myDF <- myMethylKit[myMethylKit$position %in% paste(myDMS$chr, myDMS$start, myDMS$end),]

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
gsg = goodSamplesGenes(datMeth, verbose = 3, 
                       minFraction = 3/4)#increase minFraction because too many Nas at the following steps
gsg$allOK 

# If the last statement returns TRUE, all genes have passed the cuts. If not, 
# we remove the offending genes and samples from the data
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datMeth)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datMeth)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datMeth = datMeth[gsg$goodSamples, gsg$goodGenes]
}

## Cluster the samples (in contrast to clustering genes that will come later) to see if there are any obvious
## outliers:
sampleTree = hclust(dist(datMeth), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2) ## All good!

######################################
## TRAIT DATA

## Form a data frame analogous to expression data that will hold the individuals traits
nrow(myMetaData)

## Resistance: load of parasites: myMetaData$No.Worms
## Tolerance: BCI TBC
myMetaData$Family_NUM = as.numeric(gsub("Fam", "",myMetaData$Family))

traitRows = match(rownames(datMeth), myMetaData$SampleID);
datTraits = myMetaData[traitRows,];
rownames(datTraits) = myMetaData[traitRows, "ID"];
datTraits = datTraits[c("trtG1G2_NUM", "Family_NUM", "No.Worms")] ## here keep interesting CONTINUOUS traits (no "Family", "Sex", "trtG1G2")
collectGarbage();

## We now have the methylation data in the variable datMeth, and the corresponding traits in the variable
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
abline(h=0.8,col="red")
abline(h=0.9,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
abline(h=30,col="red")
abline(h=100,col="red")
## Figure 1: Analysis of network topology for various soft-thresholding powers. The left panel shows the scale-free fit
# index (y-axis) as a function of the soft-thresholding power (x-axis). The right panel displays the mean connectivity
# (degree, y-axis) as a function of the soft-thresholding power (x-axis)
## We choose the power 5 following https://www.biostars.org/p/9477991/

#########
## 2.2 One-step network construction and module detection
net = blockwiseModules(datMeth, power = 12, # 12 is default for signed network
                       TOMType = "signed", # preserved direction +/- of correlations
                       maxBlockSize = 10000, # we have about 6000 probes, this value should be higher
                       corType = "bicor", # handles better outliers than Pearson
                       maxPOutliers = 0.1,
                       replaceMissingAdjacencies = TRUE, # to deal with pairs of genes in your data
                       # that have missing values arranged such that after removing the entries
                       # that are missing in the other gene, one of the genes becomes constant
                       saveTOMs = TRUE,
                       saveTOMFileBase = "../../data/networks/parentalTOM",
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
TOM = TOMsimilarityFromExpr(datMeth, power = 12);

rownames(TOM) <- names(datMeth)
colnames(TOM) <- names(datMeth)
# this is an adjacency matrix

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