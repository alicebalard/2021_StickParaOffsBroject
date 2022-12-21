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

if (sourceSubUnite == TRUE){
  ######################################
  ## For DMR: tile the methylation data:
  ## Kostas MBE 2020: "To identify DMRs, we used the tileMethylCounts() function in MethylKit
  # v.1.5.0 with a sliding window size of 100 bases and step size of 100 bases."
  # Summarize methylation information over tiling windows with a sliding window size of 100 bases and step size of 100 bases
  
  # 0. PARENTS trt-ctrl (CpG covered in half trt group)
  tiles_G1_half = tileMethylCounts(uniteCov6_G1_woSexAndUnknowChrOVERLAP, win.size=100,step.size=100,cov.bases = 10)
  nrow(tiles_G1_half) # methylBase object with 20348 rows
  
  # 1. CC-TC = CONTROL fish (parent CvsT) (CpG covered in half trt group)
  uniteCov14_G1bothTrt_G2control_woSexAndUnknowChr <- reorganize(methylObj = uniteCov14_G2_woSexAndUnknowChrOVERLAP,
                                                                 treatment = fullMetadata_OFFS$trtG1G2_NUM[fullMetadata_OFFS$trtG1G2_NUM %in% c(2,5)],
                                                                 sample.ids = fullMetadata_OFFS$ID[fullMetadata_OFFS$trtG1G2_NUM %in% c(2,5)])
  tiles_G1bothTrt_G2control_half = tileMethylCounts(uniteCov14_G1bothTrt_G2control_woSexAndUnknowChr,win.size=100,step.size=100,cov.bases = 10)
  nrow(tiles_G1bothTrt_G2control_half) # methylBase object with 20348 rows
  
  # 2. CT-TT = TREATMENT fish (parent CvsT) (CpG covered in half trt group)
  uniteCov14_G1bothTrt_G2infected_woSexAndUnknowChr <- reorganize(methylObj = uniteCov14_G2_woSexAndUnknowChrOVERLAP,
                                                                  treatment = fullMetadata_OFFS$trtG1G2_NUM[fullMetadata_OFFS$trtG1G2_NUM %in% c(3,6)],
                                                                  sample.ids = fullMetadata_OFFS$ID[fullMetadata_OFFS$trtG1G2_NUM %in% c(3,6)])
  tiles_G1bothTrt_G2infected_half = tileMethylCounts(uniteCov14_G1bothTrt_G2infected_woSexAndUnknowChr,win.size=100,step.size=100,cov.bases = 10)
  nrow(tiles_G1bothTrt_G2infected_half) # methylBase object with 20348 rows
  
  # 3. CC-CT = fish from CONTROL parents (G2 CvsT) (CpG covered in half trt group)
  uniteCov14_G2_woSexAndUnknowChr_G1CONTROL <- reorganize(methylObj = uniteCov14_G2_woSexAndUnknowChrOVERLAP,
                                                          treatment = fullMetadata_OFFS$trtG1G2_NUM[fullMetadata_OFFS$trtG1G2_NUM %in% c(5,6)],
                                                          sample.ids = fullMetadata_OFFS$ID[fullMetadata_OFFS$trtG1G2_NUM %in% c(5,6)])
  tiles_G2_G1CONTROL_half = tileMethylCounts(uniteCov14_G2_woSexAndUnknowChr_G1CONTROL,win.size=100,step.size=100,cov.bases = 10)
  nrow(tiles_G2_G1CONTROL_half) # methylBase object with 20348 rows
  
  # 4. TC-TT = fish from TREATMENT parents (G2 CvsT) (CpG covered in half trt group)
  uniteCov14_G2_woSexAndUnknowChr_G1INFECTED <- reorganize(methylObj = uniteCov14_G2_woSexAndUnknowChrOVERLAP,
                                                           treatment = fullMetadata_OFFS$trtG1G2_NUM[fullMetadata_OFFS$trtG1G2_NUM %in% c(2,3)],
                                                           sample.ids = fullMetadata_OFFS$ID[fullMetadata_OFFS$trtG1G2_NUM %in% c(2,3)])
  tiles_G2_G1INFECTED_half = tileMethylCounts(uniteCov14_G2_woSexAndUnknowChr_G1INFECTED,win.size=100,step.size=100,cov.bases = 10)
  nrow(tiles_G2_G1INFECTED_half) # methylBase object with 20348 rows
}
