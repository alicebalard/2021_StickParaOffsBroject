## load all libraries needed for the project

library(cowplot)
library(dendextend) # The package dendextend contains many functions for changing the appearance of a dendrogram and for comparing dendrograms.
library(devtools)
library(dplyr)
library(emmeans) ## for post-hoc Tukey tests
library(genomation) ## for annotation
# BiocManager::install("GenomicFeatures")
library(GenomicFeatures) ## for annotation
library(ggeffects) # to plot random effects predictions
library(ggplot2)
library(ggpubr) ## to merge ggplot2 plots
theme_set(theme_pubr())
library(goeveg) # find the best number of dimensions for NMDS
library(ggsignif) ## for significance bars on ggplot
library(grid)
library(gridExtra)
library(lme4) ## for mixed models
library(lmtest) # for lrtests
library(magrittr)      # provides the %>% operator
library(methylKit)
library(nlme) ## for mixed models
# install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)
library(plyr) # for join (keep row order)
library(png)
library(qualpalr)# extra palettes
library(RColorBrewer) # for colors in Venn diagrams
library(reshape2)
library(tidyverse)     # tidyverse will pull in ggplot2, readr, other useful libraries
library(vegan) ## for Adonis
library(VennDiagram)
#BiocManager::install("WGCNA")
library(WGCNA)

## offspring colors for all kind of plots
colOffs <- c("#ffe67f", "#ff6300","#a8caff","#a800d4")



