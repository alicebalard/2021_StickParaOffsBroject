## load all libraries needed for the project

list.of.packages <- c(
  "AnnotationDbi", # gene annotation from online databases
  "ape", #for reag.gff
  "biomaRt", # to retrieve genes descriptions
  "cAIC4",
  "ComplexUpset", # for prettier upset plots
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
  "ggVennDiagram",## Venn diagram in ggplot
  "goEnrichment",
  "grid",
  "gridExtra",
  "lme4", ## for mixed models
  "lmtest", # for lrtests
  "lmerTest", # for stepwise analysis of lmer
  "magrittr",      # provides the %>% operator
  "methylKit",
  "nlme", ## for mixed models
  "org.Hs.eg.db", # gene annotation from online databases
  "pairwiseAdonis",
  "pheatmap", # for heatmaps
  "plyr", # for join (keep row order",
  "png",
  "qualpalr",# extra palettes
  "RColorBrewer", # for colors in Venn diagrams
  "rentrez", # to extract info from NCBI Entrez
  "reshape2",
  "sjPlot", # plot interaction effects
  "slider", # for slidding windows
  "splitstackshape", # to spread the V9 column of gff into columns by key
  "stringr", # to modify characters
  "tidyverse",  # tidyverse will pull in ggplot2, readr, other useful libraries
  "UpSetR", # for upset plots
  "VCA",
  "vegan", ## for Adonis
  "VennDiagram")

##########################################
## install packages from github if not yet
install_github("ropensci/rentrez")
install_github("asishallab/goEnrichment")
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
install_github("gaospecial/ggVennDiagram")

###################################################################
## install from CRAN and require all libraries from CRAN and github
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
ipak(list.of.packages)

#####################################################
## install from biocmanager and require all libraries
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
