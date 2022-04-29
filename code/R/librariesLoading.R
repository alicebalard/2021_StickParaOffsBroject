## load all libraries needed for the project

list.of.packages <- c(
    "cowplot",
    "dendextend", # The package dendextend contains many functions for changing the appearance of a dendrogram and for comparing dendrograms.
    "devtools",
    "dplyr",
    "emmeans", ## for post-hoc Tukey tests
    "genomation", ## for annotation
    "ggeffects", # to plot random effects predictions
    "ggplot2",
    "ggpubr", ## to merge ggplot2 plots
    "goeveg", # find the best number of dimensions for NMDS
    "ggsignif", ## for significance bars on ggplot
    "grid",
    "gridExtra",
    "lme4", ## for mixed models
    "lmtest", # for lrtests
    "magrittr",      # provides the %>% operator
    "methylKit",
    "nlme", ## for mixed models
    "plyr", # for join (keep row order",
    "png",
    "qualpalr",# extra palettes
    "RColorBrewer", # for colors in Venn diagrams
    "reshape2",
    "splitstackshape", # to spread the V9 column of gff into columns by key
    "tidyverse",  # tidyverse will pull in ggplot2, readr, other useful libraries
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
list.bioc <- c("WGCNA",
               "GenomicFeatures") ## for annotation
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


