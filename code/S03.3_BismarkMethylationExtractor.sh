#!/bin/bash
# produces tab delimited files with content-specific methylation information
# Bismark Methylation Extractor
#$ -N bismark_methylation_extractor
#$ -o /data/scratch/btx915/04Bismark_MethylationExtraction/run_bismark_methylation_extractor.stdout
#$ -e /data/scratch/btx915/04Bismark_MethylationExtraction/run_bismark_methylation_extractor.stderr
#$ -cwd
#$ -V
#$ -l h_rt=240:00:00
#$ -l h_vmem=36G

module load bismark
module load samtools

DATA_DIR=/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/03Bismark_alignment
OUTDIR=/data/scratch/btx915/04Bismark_MethylationExtraction

cd $DATA_DIR

for file in `ls *trimmed_cutadapt_bismark_bt2.sam.gz`; 
do bismark_methylation_extractor -o $OUTDIR -s --bedGraph --counts --report $file ; done