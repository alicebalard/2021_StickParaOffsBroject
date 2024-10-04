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
  "RhpcBLASctl", # To deal with alarming nodes issue with PQLseq
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
  "VennDiagram",
  "RcppArmadillo", "foreach","parallel", "PQLseq") # all needed for PQLseq

###################################################################
message("Install CRAN packages if missing, and load CRAN packages...")

## install from CRAN and require all libraries from CRAN and github
install_if_missing <- function(packages, dependencies = TRUE) {
    new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
    if(length(new_packages)) {
        message("Installing missing packages: ", paste(new_packages, collapse = ", "))
        install.packages(new_packages, dependencies = dependencies)
    } else {
        message("All CRAN packages are already installed.\n")
    }
}

install_if_missing(list.of.packages)

message("if db error, downgrade dplyr")
# devtools::install_version("dbplyr", version = "2.3.4")

message("Loading CRAN packages...")
load_packages <- function(package_list) {
    not_loaded <- character()
    for (pkg in package_list) {
        if (suppressPackageStartupMessages(require(pkg, character.only = TRUE, quietly = TRUE))) {
                                        # Package loaded successfully, do nothing
        } else {
            not_loaded <- c(not_loaded, pkg)
        }
    }  
    if (length(not_loaded) > 0) {
        message("The following packages could not be loaded: ", paste(not_loaded, collapse = ", "))
    } else {
        message("All packages loaded successfully.\n")
    }
}

load_packages(list.of.packages)

##########################################
## install packages from github if not yet
packages_to_install <- c(
    "ropensci/rentrez","asishallab/goEnrichment","pmartinezarbizu/pairwiseAdonis/pairwiseAdonis",
    "gaospecial/ggVennDiagram", "ebbertd/chisq.posthoc.test")

install_and_load_github_packages <- function(packages) {
  # Ensure devtools is installed and loaded
  if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")
  }
  library(devtools)
  
  not_loaded <- character()
  
  for (pkg in packages) {
    pkg_name <- strsplit(pkg, "/")[[1]][2]
    
    if (!requireNamespace(pkg_name, quietly = TRUE)) {
      message(paste("Installing", pkg, "from GitHub..."))
      tryCatch(
        devtools::install_github(pkg, quiet = TRUE),
        error = function(e) {
          message(paste("Error installing", pkg, ":", e$message))
        }
      )
    }
    
    if (suppressPackageStartupMessages(require(pkg_name, character.only = TRUE, quietly = TRUE))) {
      # Package loaded successfully, do nothing
    } else {
      not_loaded <- c(not_loaded, pkg_name)
    }
  }
  
  if (length(not_loaded) > 0) {
    message("The following packages could not be loaded: ", paste(not_loaded, collapse = ", "))
  } else {
    message("All packages loaded successfully.\n")
  }
}

message("Install github packages if missing, and load github packages...")
install_and_load_github_packages(packages_to_install)

#####################################################
## install from biocmanager and require all libraries
## Biocmanager packages 
bioc_packages <- c("Category", # for hypergeometric GO test
               "WGCNA", # for networks
               "BiocFileCache", # should solve some db issues
               "ComplexHeatmap", # for heatmaps
               "genomation", ## for annotation
               "GenomicFeatures",## for annotation
               "GOstats", # for GO analysis
               "GSEABase",  # for GO term GeneSetCollection
               "GO.db", # for GO slim retrieving
               "methylKit",
               "qvalue", # for FDR after PQLseq
               "org.Hs.eg.db", # gene annotation from online databases
               "UniProt.ws" # to retrieve uniprot information
) 
# devtools::install_version("dbplyr", version = "2.3.4") # if buggy UniProt.ws

install_and_load_bioc_packages <- function(package_list) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  
  not_loaded <- character()
  
  for (pkg in package_list) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(paste("Installing", pkg, "..."))
      BiocManager::install(pkg, update = FALSE, ask = FALSE)
    }
    
    if (suppressPackageStartupMessages(require(pkg, character.only = TRUE, quietly = TRUE))) {
      # Package loaded successfully, do nothing
    } else {
      not_loaded <- c(not_loaded, pkg)
    }
  }
  
  if (length(not_loaded) > 0) {
    message("The following packages could not be loaded: ", paste(not_loaded, collapse = ", "))
  } else {
    message("All packages loaded successfully.\n")
  }
}

message("Loading bioconductor packages...")
install_and_load_bioc_packages(bioc_packages)

############# Extra configuration
## offspring colors for all kind of plots
## colOffs <- c("#ffe67f", "#ff6300","#a8caff","#a800d4") ## previous ugly
colOffs <- c("#a9d6c1ff", "#d1b000ff","#4b8da3ff","#c05a29ff")
names(colOffs) = c("NE_control", "NE_exposed", "E_control", "E_exposed")

colEffects <- as.vector(palette.colors(palette = "Okabe-Ito")[1:4])
names(colEffects) <- c("infection induced", "intergenerational", "additive", "interaction")

theme_set(theme_pubr())

########################################
## Load functions needed in the analyses
message("Load functions needed in the analyses...")
source("customRfunctions.R")
source("homebrewDMSannotation.R")
message("done.\n")

###################
## Annotation files
message("Load annotation files: GYgynogff containing length of each gynogen chromosomes, annotBed12 and annotGff3 for the full annotation...")

## Load file containing length of each gynogen chromosomes
## grep "contig" gitignore/bigdata/Gy_allnoM_rd3.maker_apocrita.noseq_corrected.gff | awk '{print $1, $5}' > data/Gy_allnoM_rd3.maker_apocrita.noseq_corrected_chromoAndLength.txt
GYgynogff = read.table("../../data/Gy_allnoM_rd3.maker_apocrita.noseq_corrected_chromoAndLength.txt")
names(GYgynogff) = c("chrom","length")
#GA_genome.fa.sizes.txt is a file with chromosome sizes and names
GYgynogff <- GYgynogff[!GYgynogff$chrom %in% c("Gy_chrUn", "M"),] %>% ## rm unknown and mitochondrial chromosome
  mutate(chrom_nr=chrom %>% deroman(),
         chrom_order=factor(chrom_nr) %>% as.numeric()) %>% 
  arrange(chrom_order) %>% # keep chr in order
  mutate(gstart=lag(length,default=0) %>% cumsum(),
         gend=gstart+length, 
         type=LETTERS[2-(chrom_order%%2)],
         gmid=(gstart+gend)/2,
         chrom=as.factor(chrom),
         type=rep(c("A","B"),length(length)/2))

## NB Promoters are defined by options at genomation::readTranscriptFeatures function. 
## The default option is to take -1000,+1000bp around the TSS and you can change that. 
## -> following Heckwolf 2020 and Sagonas 2020, we consider 1500bp upstream and 500 bp downstream
annotBed12=readTranscriptFeatures(
  "../../gitignore/bigdata/06GynoAnnot/Gy_allnoM_rd3.maker_apocrita.noseq_corrected.gff.streamlined_for_AGAT.gff.fun.bed12",
  remove.unusual = FALSE, up.flank = 1500, down.flank = 500, unique.prom = F) ## unique.prom otherwise promoters have no name

## Load curated gff file
annotGff3 <- rtracklayer::readGFF("../../gitignore/bigdata/06GynoAnnot/Gy_allnoM_rd3.maker_apocrita.noseq_corrected.gff.streamlined_for_AGAT.gff.fun.gff3")

message("done.\n")
