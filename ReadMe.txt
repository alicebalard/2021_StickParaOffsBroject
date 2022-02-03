## Alice Balard
## 25th of May 2021

############## DATA ORGANISATION  ################
[btx915@frontend8:/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/]
├── Data
│   ├── 01TrimmedReads_cutadapt
│   ├── 02RefGenome
│   └── 03Bismark_alignment
├── GIT_StickParaOffsBroject
│   ├── code
│   ├── data
│   └── ReadMe.txt
└── QualityControl
    ├── QC00_RawReads_fastqc
    └── QC01_TrimmedReads_fastqc

########## SCRIPTS (within GIT_StickParaOffsBroject/code) ################
├── calculateStats.R
├── ReadMe.txt
├── S00_multiQC_generateReports.sh
├── S01_fastqc.sh
├── S02.1_adapterTrim_cutadapt.sh
├── S02.2_fastqcTrimmed.sh
├── S03.1_GenomePreparation.sh
└── S03.2_BismarkAlignment

############# STEPS ##############

### Raw files fro the 114 samples (zipped) here:
/data/archive/archive-SBCS-EizaguirreLab/RRBS/StickPara_Broject/00Illumina_RawReads/

# Quality control of the raw sequence reads from the bisulfite-treated samples was done with FASTQC v0.11.9 (Andrews 2010) 
# (and after trimming to check improvement)
Code: S01_fastqc.sh
Output: QualityControl/QC00_RawReads_fastqc 
        QualityControl/QC01_TrimmedReads_fastqc

# Run multiQC v1.10.1 for having one report (will be done at different steps hence the 00)
Code: S00_multiQC_generateReports.sh
Output: QC00_RawReads_fastqc/multiqc_report.html
        QC01_TrimmedReads_fastqc/multiqc_report.html

# Reads were processed and filtered to remove adaptor sequences and low-quality (i.e., q < 20) reads with Cutadapt v1.13 (Martin 2011) 
# using three adapter sequences (NNAGATCGGAAGAGCACAC, AGATCGGAAGAGCACAC, ATCGGAAGAGCACAC).
Code: S02.1_adapterTrim_cutadapt.sh
Output: 01TrimmedReads_cutadapt

# Reference gynogen genome preparation
Code: S03.1_GenomePreparation.sh
Output: 02RefGenome/Bisulfite_Genome

# Bismark alignment 
Code: S03.2_BismarkAlignment
Output: 03Bismark_alignment

##############################
##### R methylation analysis workflow: #####
Preparation of data and metadata:
R00 -> alignment stats
R01.1, R01.2, R01.3, R01.4 -> data preparation

R02 -> link methylation and fitness (tbc)

# PART 1. (script R03.1) Methylation profile: or Adonis multivariate stats, presence/abs, block effect=family
# PART 2. (script R03.2) DMS 
# PART 3.DMR: Precise localisation of blocks investing in methylation
# PART 4.Network: Explore other parameters than geographic: are the modules
#   on a high recombination place? Or with high mutation rate? 
#   What does modules capture on top of DMRs?
##### Then 
# - Gene Ontology
# - Final: Check robustness: reshuffle the controls and check that we find less DMS/DMR/modules
# 
# These analyses will be done on:
# A. Parents ctr-trt
# B. Offspring ctr-trt pat1 & ctrl-trt pat2 (different offspring treatments)
# C. Offspring ctrl-ctrl & trt-trt (different paternal treatment)
##############################
