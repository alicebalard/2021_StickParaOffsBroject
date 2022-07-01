## load libraries
source("librariesLoading.R")
## load custom functions
source("customRfunctions.R")
## Load samples metadata
source("R02.1_loadMetadata.R")
## Load previously united methylkit data
source("R02.2_loadMethyldata.R") 

if (loadannot == TRUE){
  #########################
  ## Load file containing length of each gynogen chromosomes
  ## grep "contig" gitignore/bigdata/Gy_allnoM_rd3.maker_apocrita.noseq_corrected.gff | awk '{print $1, $5}' > data/Gy_allnoM_rd3.maker_apocrita.noseq_corrected_chromoAndLength.txt
  GYgynogff = read.table("../../data/Gy_allnoM_rd3.maker_apocrita.noseq_corrected_chromoAndLength.txt")
  names(GYgynogff) = c("chrom","length")
  #########################
  ## Load genome annotation
  ## NB Promoters are defined by options at genomation::readTranscriptFeatures function. 
  ## The default option is to take -1000,+1000bp around the TSS and you can change that. 
  ## -> following Heckwolf 2020 and Sagonas 2020, we consider 1500bp upstream and 500 bp downstream
  annotBed12=readTranscriptFeatures("../../gitignore/bigdata/06GynoAnnot/Gy_allnoM_rd3.maker_apocrita.noseq_corrected.gff.streamlined_for_AGAT.CURATED.transdec.bed12",
                                    remove.unusual = FALSE, up.flank = 1500, down.flank = 500)
  
  ## Change recursively the gene names to keep only ID
  getName <- function(x) {sub(";.*", "", sub(".*ID=", "", x))}
  for (i in 1:length(annotBed12)){
    annotBed12[[i]]$name <- getName(annotBed12[[i]]$name)
  }
  
  ## Load curated gff file
  annotGff3 <- rtracklayer::readGFF("../../gitignore/bigdata/06GynoAnnot/Gy_allnoM_rd3.maker_apocrita.noseq_corrected.gff.streamlined_for_AGAT.CURATED.gff")
}

## That can only be done AFTER these have been calculated of course. Only for later scripts.
if (sourceDMS == TRUE){
  ## Source the previously calculated DMS/DMR
  DMSlist <- readRDS("../../data/DiffMeth/DMSlist.RDS")
  DMRlist <- readRDS("../../data/DiffMeth/DMRlist.RDS")
  ## By BP
  DMBPlist <- readRDS("../../data/DiffMeth/DMperBP_list.RDS")
  DMSBPlist <- lapply(DMBPlist, "[[", 1)
  DMRBPlist <- lapply(DMBPlist, "[[", 2)
}
