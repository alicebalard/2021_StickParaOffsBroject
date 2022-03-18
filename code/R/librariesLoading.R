## load all libraries needed for the project
library(plyr) # for join (keep row order)
library(dendextend) # The package dendextend contains many functions for changing the appearance of a dendrogram and for comparing dendrograms.
library(ggplot2)
library(cowplot)
library(ggsignif) ## for significance bars on ggplot
library(lme4) ## for mixed models
library(nlme) ## for mixed models
library(tidyverse)
library(emmeans) ## for post-hoc Tukey tests
library(methylKit)
library(vegan) ## for Adonis
library(genomation) ## for annotation
# BiocManager::install("GenomicFeatures")
library(GenomicFeatures) ## for annotation
library(tidyverse)     # tidyverse will pull in ggplot2, readr, other useful libraries
library(magrittr)      # provides the %>% operator
#BiocManager::install("WGCNA")
library(WGCNA)
library(ggpubr) ## to merge ggplot2 plots
theme_set(theme_pubr())
library(goeveg) # find the best number of dimensions for NMDS
library(VennDiagram)
library(RColorBrewer) # for colors in Venn diagrams
library(lmtest) # for lrtests
library(ggeffects) # to plot random effects predictions
library(reshape2)
library(png)
library(grid)

## offspring colors for all kind of plots
colOffs <- c("#ffe67f", "#ff6300","#a8caff","#a800d4")



