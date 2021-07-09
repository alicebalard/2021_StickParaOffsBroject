#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=1G
#$ -l h_rt=240:0:0
#$ -cwd
#$ -j y

#Run fastqc for quality check of raw reads

module load fastqc

cd /data/scratch/btx915/00Illumina_RawReads/

fastqc *fastq.gz

# move files to the result directory after run
mv *fastqc* /data/SBCS-EizaguirreLab/Alice/StickParaBroOff/QualityControl/QC00_RawReads_fastqc

