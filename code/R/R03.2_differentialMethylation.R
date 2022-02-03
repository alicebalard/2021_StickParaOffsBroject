## Differential methylation analyses
## A. Balard
## February 2022

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

###############################################
### PART 2: Differential methylation sites ####
###############################################

##########################################################
## Comparison 1: BASELINE -> Parents (control vs infected) 
## Extract DMS from parents (at least in 2 fish), annotate them, compare with Kostas results

## Calculate DMS accounting for covariates: family
cov = data.frame(Family = fullMetadata_PAR$Family)
myDiffMeth=calculateDiffMeth(uniteCov2_woSexAndUnknowChr_PAR, 
                             covariates = cov, mc.cores = 4)

# We select the bases that have q-value<0.01 and percent methylation difference larger than 15%.
# NB: arg type="hyper" or type="hypo" gives hyper-methylated or hypo-methylated regions/bases.
myDiff1_15p = getMethylDiff(myDiffMeth,difference=15,qvalue=0.01)

myDiff1_15p # 6544 positions
saveRDS(myDiff1_15p, file = "../../gitignore/output/myDiff1_15p_parentalDiffMeth.RDS")
# 
# myDiff1_15p <- readRDS("../../gitignore/output/myDiff1_15p_parentalDiffMeth.RDS")
# 
# # annotation!!
# library(genomation)
# 
# #library(GenomicFeatures)
# 
# gene.obj=readTranscriptFeatures("../../gitignore/bigdata/Gy_allnoM_rd3.maker_apocrita.noseq_corrected.bed12", remove.unusual = FALSE)
# 
# # annotate differentially methylated CpGs with promoter/exon/intron using annotation data
# annotateWithGeneParts(as(myDiff1_15p,"GRanges"),gene.obj)
# 
# ## Kostas MBE: The DMSs and regions were predominately
# #found in intergenic regions (47.74% and 48.94%, respecti vely),with introns (26.19% and 23.09), exons (15.07% and 13.98%),and promoters (11% and 13.98%) showing lower proportions
# 
# ?annotateWithGeneParts()
# 
# 
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
# #### Annotating differentially methylated bases or regions
# #BiocManager::install("genomation")
# # library(genomation)
# ## TBC
# 
# ##########################################################
# #### Co-methylation network construction
# ### Help: https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-01-dataInput.pdf
# library(tidyverse)     # tidyverse will pull in ggplot2, readr, other useful libraries
# library(magrittr)      # provides the %>% operator
# #BiocManager::install("WGCNA") 
# library(WGCNA) 
# 
# #######
# ## 1.Data preparation
# 
# ## Extract the targeted CpG positions from the methylkit object
# myMethylKit$position <- paste(myMethylKit$chr, myMethylKit$start, myMethylKit$end)
# myDF <- myMethylKit[myMethylKit$position %in% paste(myDiff1_15p$chr, myDiff1_15p$start, myDiff1_15p$end),]
# 
# ## Make myDF_fr : a dataframe in which each column is a sample or auxilliary info,
# ## and each row a CpG position, with values = beta value (methylation frequency = numCs/coverage)
# myDF <-methylKit::getData(myDF)
# myDF_fr <- data.frame(matrix(ncol=length(myMethylKit@sample.ids), nrow=nrow(myDF)))
# names(myDF_fr) = myMethylKit@sample.ids
# for(i in 1:ncol(myDF_fr)){
#   myDF_fr[i] <- myDF[paste0("numCs", i)]/myDF[paste0("coverage", i)]
# }
# 
# ## Store CpG information (Chr, start, end)
# CpGInfo <- data.frame(CpGpos = paste0("CpG", 1:length(myDF$chr)), 
#                       Chromosome = myDF$chr,
#                       start = myDF$start,
#                       end=myDF$end)
# 
# ## We now transpose the expression data for further analysis.
# datMeth = as.data.frame(t(myDF_fr));
# names(datMeth) = CpGInfo$CpGpos
# 
# ## We first check for genes and samples with too many missing values:
# gsg = goodSamplesGenes(datMeth, verbose = 3);
# gsg$allOK # all good!
# 
# ## Cluster the samples (in contrast to clustering genes that will come later) to see if there are any obvious
# ## outliers:
# sampleTree = hclust(dist(datMeth), method = "average");
# # Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# # The user should change the dimensions if the window is too large or too small.
# sizeGrWindow(12,9)
# #pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
# par(cex = 0.6);
# par(mar = c(0,4,2,0))
# plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
#      cex.axis = 1.5, cex.main = 2) ## All good!
# 
# ## Form a data frame analogous to expression data that will hold the individuals traits
# 
# nrow(myMetaData)
# 
# traitRows = match(rownames(datMeth), myMetaData$SampleID);
# datTraits = myMetaData[traitRows,];
# rownames(datTraits) = myMetaData[traitRows, "ID"];
# datTraits = datTraits[c("BCI", "integration.RLU.")] ## here keep interesting CONTINUOUS traits (no "Family", "Sex", "trtG1G2")
# collectGarbage();
# 
# ## We now have the expression data in the variable datMeth, and the corresponding traits in the variable
# ## datTraits. Before we continue with network construction and module detection, we visualize how 
# ## the traits relate to the sample dendrogram
# # Re-cluster samples
# sampleTree2 = hclust(dist(datMeth), method = "average")
# # Convert traits to a color representation: white means low, red means high, grey means missing entry
# traitColors = numbers2colors(datTraits, signed = FALSE);
# # Plot the sample dendrogram and the colors underneath.
# plotDendroAndColors(sampleTree2, traitColors,
#                     groupLabels = names(datTraits),
#                     main = "Sample dendrogram and trait heatmap")
# 
# ## if needed
# ## save(datMeth, datTraits, file = "filename.RData")
# 
# #######
# ## 2.Automatic network construction and module detection
# 
# ###
# ## 2.1 Choosing the soft-thresholding power: analysis of network topology
# # Choose a set of soft-thresholding powers
# powers = c(c(1:10), seq(from = 12, to=20, by=2))
# # Call the network topology analysis function
# sft = pickSoftThreshold(datMeth, powerVector = powers, verbose = 5)
# # Plot the results:
# sizeGrWindow(9, 5) ; par(mfrow = c(1,2)); cex1 = 0.9;
# # Scale-free topology fit index as a function of the soft-thresholding power
# plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
#      xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
#      main = paste("Scale independence"));
# text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
#      labels=powers,cex=cex1,col="red");
# # this line corresponds to using an R^2 cut-off of h
# abline(h=0.90,col="red")
# # Mean connectivity as a function of the soft-thresholding power
# plot(sft$fitIndices[,1], sft$fitIndices[,5],
#      xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
#      main = paste("Mean connectivity"))
# text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
# ## Figure 1: Analysis of network topology for various soft-thresholding powers. The left panel shows the scale-free fit
# # index (y-axis) as a function of the soft-thresholding power (x-axis). The right panel displays the mean connectivity
# # (degree, y-axis) as a function of the soft-thresholding power (x-axis)
# ## We choose the power 3, which is the lowest power for which the scale-free topology fit
# # index curve flattens out upon reaching a high value.
# 
# ###
# ## 2.2 One-step network construction and module detection
# net = blockwiseModules(datMeth, power = 3,
#                        TOMType = "unsigned", minModuleSize = 30,
#                        reassignThreshold = 0, mergeCutHeight = 0.25,
#                        numericLabels = TRUE, pamRespectsDendro = FALSE,
#                        saveTOMs = TRUE,
#                        saveTOMFileBase = "femaleMouseTOM",
#                        verbose = 3)
# 
# # We now return to the network analysis. To see how many modules were identified and what the module sizes are,
# # one can use table(net$colors). Its output is
# table(net$colors) # 1 module, with 42 CpG. The label 0 is reserved for CpGs outside of all modules
# 
# # open a graphics window
# sizeGrWindow(12, 9)
# # Convert labels to colors for plotting
# mergedColors = labels2colors(net$colors)
# # Plot the dendrogram and the module colors underneath
# plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
#                     "Module colors",
#                     dendroLabels = FALSE, hang = 0.03,
#                     addGuide = TRUE, guideHang = 0.05)
# 
# moduleLabels = net$colors
# moduleColors = labels2colors(net$colors)
# MEs = net$MEs;
# geneTree = net$dendrograms[[1]];
# # save(MEs, moduleLabels, moduleColors, geneTree,
# #      file = "FemaleLiver-02-networkConstruction-auto.RData")
# 
# ################ GRAPHS ###################
# library(tidyverse); library(igraph); library(ggraph); library(tidygraph); library(ggregplot)
# library(fs); library(magrittr); library(cowplot); library(ggplot2)
# library(patchwork); library(Matrix); library(MASS); library(reshape2)
# 
# ## Export from WGCNA to igraph for visualisation
# 
# # Recalculate topological overlap
# TOM = TOMsimilarityFromExpr(datMeth, power = 3);
# 
# rownames(TOM) <- names(datMeth)
# colnames(TOM) <- names(datMeth)
# ## 66 by 66: this is an adjacency matrix
# 
# ## Make graph:
# myGraph <- TOM %>% 
#   # stores the correlations as edge weight:
#   graph_from_adjacency_matrix(mode = "undirected", weighted = T) %>% # igraph. Weight saved in edges
#   # as_tbl_graph(nodes = CpGInfo) %>% # %>% # tidygraph with nodes CpG info
#   simplify()%>% # removes self-loops
#   as_tbl_graph() %>% # tidygraph with nodes CpG name
#   activate(nodes) %>% 
#   left_join(CpGInfo, by = c("name" = "CpGpos"))%>% # add CpGInfo: chromosome location 
#   left_join(data.frame(name= names(moduleLabels), color = moduleColors), by = c("name" = "name")) # add previously detected module colors
# 
# ## add layout based on algorythms, compare different layouts
# Layout1 <- layout_with_fr(myGraph, # Fruchterman-Reingold as an example
#                           weights = myGraph %>% 
#                             activate(edges) %>% 
#                             pull(weight)) # Adding in weights explicitly
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
# 
# ## Plot:
# myGraph %>% # my igraph - tidygraph
#   ggraph(Layout1) +  # then ggraph
#   geom_edge_link(aes(alpha=weight)) + 
#   scale_edge_alpha(range=c(0.1,1))+
#   geom_node_point(aes(color=color), size = 7) +
#   scale_colour_manual(values = c("grey", "turquoise"))+
#   geom_node_text(aes(label = name), size = 3) 
# 
# 
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
# 
