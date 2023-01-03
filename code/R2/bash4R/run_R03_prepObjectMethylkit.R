#!/bin/bash
#$ -pe smp 10
#$ -l h_vmem=10G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -j y
#$ -N runR
#$ -o /data/SBCS-EizaguirreLab/Alice/StickParaBroOff/StickParaOffsBroject/code/R2/bash4R/run_R03_prepObjectMethylkit.R.stdout
#$ -e /data/SBCS-EizaguirreLab/Alice/StickParaBroOff/StickParaOffsBroject/code/R2/bash4R/run_R03_prepObjectMethylkit.R.stderr

module load R/4.0.2

cd /data/SBCS-EizaguirreLab/Alice/StickParaBroOff/StickParaOffsBroject/code/R2/

Rscript /data/SBCS-EizaguirreLab/Alice/StickParaBroOff/StickParaOffsBroject/code/R2/R03_prepObjectMethylkit_runInCLUSTER.R
