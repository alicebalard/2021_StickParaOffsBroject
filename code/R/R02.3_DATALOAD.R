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
  annotBed12=readTranscriptFeatures("../../gitignore/bigdata/06GynoAnnot/Gy_allnoM_rd3.maker_apocrita.noseq_corrected.gff.streamlined_for_AGAT.CURATED.bed12",
                                    remove.unusual = FALSE, up.flank = 1500, down.flank = 500)
  
  ## Load curated gff file
  annotGff3 <- rtracklayer::import("../../gitignore/bigdata/06GynoAnnot/Gy_allnoM_rd3.maker_apocrita.noseq_corrected.gff.streamlined_for_AGAT.CURATED.gff")
}

## That can only be done AFTER these have been calculated of course. Only for later scripts.
if (sourceDMS == TRUE){
  ## Source the previously calculated DMS/DMR
  ## Parents (brotherPairID as covariates)
  ### DM from CpG positions shared by half the fish per trt
  DMS15pc_G1_half <- readRDS("../../data/DiffMeth/DMS15pc_G1_half.RDS"); nrow(DMS15pc_G1_half) # 3648
  DMR15pc_G1_half <- readRDS("../../data/DiffMeth/DMR15pc_G1_half.RDS"); nrow(DMR15pc_G1_half) # 23
  ### DM from CpG positions shared by all the fish
  # DMS15pc_G1_ALL <- readRDS("../../data/DiffMeth/DMS15pc_G1_ALL.RDS"); nrow(DMS15pc_G1_ALL) # 125
  # DMR15pc_G1_ALL returned 0 DMR
  
  ## Offspring (brotherPairID & Sex as covariates)
  ## Control G1 - G2(trt vs control)
  ### DM from CpG positions shared by half the fish per trt
  DMS15pc_G2_controlG1_half <- readRDS("../../data/DiffMeth/DMS15pc_G2_controlG1_half.RDS")
  nrow(DMS15pc_G2_controlG1_half) # 1197
  DMR15pc_G2_controlG1_half <- readRDS("../../data/DiffMeth/DMR15pc_G2_controlG1_half.RDS")
  nrow(DMR15pc_G2_controlG1_half) # 6
  ### DM from CpG positions shared by all the fish
  # DMS15pc_G2_controlG1_ALL <- readRDS("../../data/DiffMeth/DMS15pc_G2_controlG1_ALL.RDS")
  # nrow(DMS15pc_G2_controlG1_ALL) # 38
  # DMR15pc_G2_controlG1_ALL returned 0 DMR
  
  ## Infected G1 - G2(trt vs control)
  ### DM from CpG positions shared by half the fish per trt
  DMS15pc_G2_infectedG1_half <- readRDS("../../data/DiffMeth/DMS15pc_G2_infectedG1_half.RDS")
  nrow(DMS15pc_G2_infectedG1_half) # 690
  DMR15pc_G2_infectedG1_half <- readRDS("../../data/DiffMeth/DMR15pc_G2_infectedG1_half.RDS")
  nrow(DMR15pc_G2_infectedG1_half) # 8
  ### DM from CpG positions shared by all the fish
  # DMS15pc_G2_infectedG1_ALL <- readRDS("../../data/DiffMeth/DMS15pc_G2_infectedG1_ALL.RDS")
  # nrow(DMS15pc_G2_infectedG1_ALL) # 22
  # DMR15pc_G2_infectedG1_ALL <- readRDS("../../data/DiffMeth/DMR15pc_G2_infectedG1_ALL.RDS")
  # nrow(DMR15pc_G2_infectedG1_ALL) # 1
  
  ## Both trt G1 - Control G2
  ### DM from CpG positions shared by half the fish per trt
  DMS15pc_G1_controlG2_half <- readRDS("../../data/DiffMeth/DMS15pc_G1_controlG2_half.RDS")
  nrow(DMS15pc_G1_controlG2_half) # 1569
  DMR15pc_G1_controlG2_half <- readRDS("../../data/DiffMeth/DMR15pc_G1_controlG2_half.RDS")
  nrow(DMR15pc_G1_controlG2_half) # 14
  
  ## Both trt G1 - Infected G2
  ### DM from CpG positions shared by half the fish per trt
  DMS15pc_G1_infectedG2_half <- readRDS("../../data/DiffMeth/DMS15pc_G1_infectedG2_half.RDS")
  nrow(DMS15pc_G1_infectedG2_half) # 2050
  DMR15pc_G1_infectedG2_half <- readRDS("../../data/DiffMeth/DMR15pc_G1_infectedG2_half.RDS")
  nrow(DMR15pc_G1_infectedG2_half) # 19
}