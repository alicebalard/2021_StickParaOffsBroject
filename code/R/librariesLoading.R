## load all libraries needed for the project

list.of.packages <- c(
  "ape", #for reag.gff
  "cowplot",
  "dendextend", # The package dendextend contains many functions for changing the appearance of a dendrogram and for comparing dendrograms.
  "devtools",
  "dplyr",
  "emmeans", ## for post-hoc Tukey tests
  "forcats", # keeps characters in previous order axis ggplot (for bubble plot)
  "genomation", ## for annotation
  "ggeffects", # to plot random effects predictions
  "ggplot2",
  "ggpubr", ## to merge ggplot2 plots
  "ggrepel",
  "goeveg", # find the best number of dimensions for NMDS
  "ggsignif", ## for significance bars on ggplot
  "grid",
  "gridExtra",
  "lme4", ## for mixed models
  "lmtest", # for lrtests
  "lmerTest", # for stepwise analysis of lmer
  "magrittr",      # provides the %>% operator
  "methylKit",
  "nlme", ## for mixed models
  "pheatmap", # for heatmaps
  "plyr", # for join (keep row order",
  "png",
  "qualpalr",# extra palettes
  "RColorBrewer", # for colors in Venn diagrams
  "reshape2",
  "sjPlot", # plot interaction effects
  "slider", # for slidding windows
  "splitstackshape", # to spread the V9 column of gff into columns by key
  "cAIC4",
  "tidyverse",  # tidyverse will pull in ggplot2, readr, other useful libraries
  "UpSetR", # for upset plots
  "VCA",
  "vegan", ## for Adonis
  "VennDiagram")

ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
ipak(list.of.packages)

## github packages manual check
if (!"pairwiseAdonis" %in% installed.packages()[, "Package"]) {
  devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
}
library(pairwiseAdonis)

## Biocmanager packages 
list.bioc <- c("Category", # for hypergeometric GO test
               "WGCNA", # for networks
               "GenomicFeatures",## for annotation
               "GOstats", # for GO analysis
               "GSEABase"  # for GO term GeneSetCollection
) 
ipak2 <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    BiocManager::install(new.pkg)
  sapply(pkg, require, character.only = TRUE)
}
ipak2(list.bioc)

## offspring colors for all kind of plots
colOffs <- c("#ffe67f", "#ff6300","#a8caff","#a800d4")

theme_set(theme_pubr())


