#!/bin/bash
#$ -pe smp 10
#$ -l h_vmem=10G
#$ -l h_rt=240:0:0
#$ -cwd
#$ -j y
#$ -N runR

module load R/4.0.2

Rscript /data/SBCS-EizaguirreLab/Alice/StickParaBroOff/StickParaOffsBroject/code/R/R04.1_getDifferentialMethylation.R
