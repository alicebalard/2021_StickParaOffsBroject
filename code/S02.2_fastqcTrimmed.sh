#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=1G
#$ -l h_rt=240:0:0
#$ -cwd
#$ -j y
#$ -o /data/SBCS-EizaguirreLab/Alice/StickParaBroOff/00Quality_Control_fastqc/QCtrimmedReads_cutadapt/run_fastqc_Alice.stdout 
#$ -e /data/SBCS-EizaguirreLab/Alice/StickParaBroOff/00Quality_Control_fastqc/QCtrimmedReads_cutadapt/run_fastqc_Alice.stderr

#Run fastqc for quality check of raw reads

module load fastqc

cd /data/scratch/btx915/01TrimmedReads_cutadapt/

fastqc *trimmed_cutadapt.fastq.gz
