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
  "factoextra", # color PCA plots
  "FactoMineR", # for PCA
  "forcats", # keeps characters in previous order axis ggplot (for bubble plot)
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
  "MuMIn", # participation of variables to the variance
  "missMDA",# PCA for incomplete data
  "nlme", ## for mixed models
  "pairwiseAdonis",
  "pheatmap", # for heatmaps
  "plyr", # for join (keep row order",
  "png",
  "purrr",
  "PQLseq", # calculate diff meth with relatedness
  "qualpalr",# extra palettes
  "RColorBrewer", # for colors in Venn diagrams
  "rcompanion", # for Spearman's rho 95%CI by BS
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

###################################################################
## install from CRAN and require all libraries from CRAN and github
install_if_missing <- function(packages, dependencies = TRUE) {
  new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
  if(length(new_packages)) {
    message("Installing missing packages: ", paste(new_packages, collapse = ", "))
    install.packages(new_packages, dependencies = dependencies)
  } else {
    message("All CRAN packages are already installed.")
  }
}

install_if_missing(list.of.packages)

message("Loading CRAN packages...")
for(pkg in list.of.packages) {
    if(!require(pkg, character.only = TRUE)) {
        warning(paste("Failed to load package:", pkg))
    }
}

##########################################
## install packages from github if not yet
packages_to_install <- c(
    "ropensci/rentrez","asishallab/goEnrichment","pmartinezarbizu/pairwiseAdonis/pairwiseAdonis", "gaospecial/ggVennDiagram")

install_and_load_github_packages <- function(packages, dependencies = TRUE) {
  # Ensure devtools is installed and loaded
  if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")
  }
  library(devtools)
  
  # Function to extract package name from GitHub repo string
  extract_package_name <- function(repo) {
    parts <- strsplit(repo, "/")[[1]]
    if (length(parts) > 1) {
      return(parts[2])
    } else {
      return(repo)
    }
  }
  
  # Install and load packages
  for (pkg in packages) {
    pkg_name <- extract_package_name(pkg)
    if (!requireNamespace(pkg_name, quietly = TRUE)) {
      message(paste("Installing", pkg, "from GitHub..."))
      tryCatch(
        install_github(pkg, dependencies = dependencies),
        error = function(e) {
          message(paste("Error installing", pkg, ":", e$message))
          return()  # Skip to next package if installation fails
        }
      )
    } else {
      message(paste(pkg_name, "is already installed."))
    }
    
    # Attempt to load the package
    message(paste("Loading", pkg_name, "..."))
    if (!require(pkg_name, character.only = TRUE)) {
      warning(paste("Failed to load package:", pkg_name))
    }
  }
}

message("Loading github packages...")
install_and_load_github_packages(packages_to_install)

#####################################################
## install from biocmanager and require all libraries
## Biocmanager packages 
bioc_packages <- c("Category", # for hypergeometric GO test
               "WGCNA", # for networks
               "genomation", ## for annotation
               "GenomicFeatures",## for annotation
               "GOstats", # for GO analysis
               "GSEABase",  # for GO term GeneSetCollection
               "methylKit",
               "org.Hs.eg.db" # gene annotation from online databases
) 

install_and_load_bioc_packages <- function(packages) {
  # Ensure BiocManager is installed and loaded
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  library(BiocManager)

  # Install and load packages
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(paste("Installing", pkg, "..."))
      BiocManager::install(pkg, update = FALSE, ask = FALSE)
    }
    
    message(paste("Loading", pkg, "..."))
    if (!require(pkg, character.only = TRUE)) {
      warning(paste("Failed to load package:", pkg))
    }
  }
}

message("Loading bioconductor packages...")
install_and_load_bioc_packages(bioc_packages)

############# Extra configuration
## offspring colors for all kind of plots
colOffs <- c("#ffe67f", "#ff6300","#a8caff","#a800d4")

theme_set(theme_pubr())
