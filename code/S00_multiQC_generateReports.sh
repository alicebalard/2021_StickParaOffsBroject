#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=1G
#$ -l h_rt=240:0:0
#$ -cwd
#$ -j y

# run from conda environment "myenv" where multiqc is loaded (conda activate myenv)

# Run multiqc on raw reads fastqc quality checks
multiqc /data/SBCS-EizaguirreLab/Alice/StickParaBroOff/01Quality_Control_fastqc/
