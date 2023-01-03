#!/bin/bash
#$ -pe smp 10
#$ -l h_vmem=10G
#$ -l h_rt=240:0:0
#$ -cwd
#$ -j y
#$ -N runR
#$ -o /data/SBCS-EizaguirreLab/Alice/StickParaBroOff/StickParaOffsBroject/code/R2/bash4R/run_R04_getDifferentialMethylation_runInCLUSTER.R.stdout
#$ -e /data/SBCS-EizaguirreLab/Alice/StickParaBroOff/StickParaOffsBroject/code/R2/bash4R/run_R04_getDifferentialMethylation_runInCLUSTER.R.stderr

module load R/4.2.0

cd /data/SBCS-EizaguirreLab/Alice/StickParaBroOff/StickParaOffsBroject/code/R2/

Rscript /data/SBCS-EizaguirreLab/Alice/StickParaBroOff/StickParaOffsBroject/code/R2/R04_getDifferentialMethylation_runInCLUSTER.R

