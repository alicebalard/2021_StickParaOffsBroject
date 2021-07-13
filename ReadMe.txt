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
