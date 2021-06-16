#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=1G
#$ -l h_rt=240:0:0
#$ -cwd
#$ -j y

#Run fastqc for quality check of raw reads

module load fastqc

cd /data/scratch/btx915/01Trimmed_Reads/

fastqc *trimmed.fq.gz

# move files to the result directory after run
mv *fastqc* /data/SBCS-EizaguirreLab/Alice/StickParaBroOff/00Quality_Control_fastqc/QCtrimmedReads/.
