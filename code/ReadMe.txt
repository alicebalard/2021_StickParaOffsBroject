## Alice Balard
## 25th of May 2021
## From Sagonas et al. 2020

### Raw files fro the 114 samples (zipped) here:
/data/archive/archive-SBCS-EizaguirreLab/RRBS/StickPara_Broject/00Illumina_RawReads/

# Quality control of the raw sequence reads from the bisulfite-treated samples was done with FASTQC v0.11.9 (Andrews 2010)
Code: S01_fastqc.sh
Output: /data/SBCS-EizaguirreLab/Alice/StickParaBroOff/01Quality_Control_fastqc/

# Run multiQC v1.10.1 for having one report (will be done at different steps)
Code: multiqc .
Output: /data/SBCS-EizaguirreLab/Alice/StickParaBroOff/01Quality_Control_fastqc/multiqc_report.html

# Reads were processed and filtered to remove adaptor sequences and low-quality (i.e., q < 20) reads with Cutadapt v1.13 (Martin 2011) using three adapter sequences (NNAGATCGGAAGAGCACAC, AGATCGGAAGAGCACAC, ATCGGAAGAGCACAC).

trimgalore/0.6.5(default) on apocrita
cutadapt (conda loaded) v1.18





# in 02Trimmed_Reads/Cutadapt cut in A, B and C folder 18 + 18 + 17 = 53 samples
