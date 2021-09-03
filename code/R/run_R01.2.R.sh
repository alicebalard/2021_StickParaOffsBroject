#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=150G
#$ -l highmem
#$ -l h_rt=1:0:0
#$ -cwd
#$ -j y
#$ -N runR
#$ -o /data/SBCS-EizaguirreLab/Alice/StickParaBroOff/GIT_StickParaOffsBroject/code/R/run_R01.2.R.stdout
#$ -e /data/SBCS-EizaguirreLab/Alice/StickParaBroOff/GIT_StickParaOffsBroject/code/R/run_R01.2.R.stderr

module load R/4.0.2

Rscript /data/SBCS-EizaguirreLab/Alice/StickParaBroOff/GIT_StickParaOffsBroject/code/R/R01.2_prepObjectMethylkit.R
