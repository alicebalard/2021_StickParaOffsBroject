#!/bin/bash
#$ -pe smp 10
#$ -l h_vmem=5G
#$ -l h_rt=240:0:0
#$ -cwd
#$ -j y
#$ -N runR04
#$ -m e
#$ -o /data/SBCS-EizaguirreLab/Alice/StickParaBroOff/StickParaOffsBroject/code/R/bash4R/run_R04_getDiffMeth_PQLseq_G1.stdout
#$ -e /data/SBCS-EizaguirreLab/Alice/StickParaBroOff/StickParaOffsBroject/code/R/bash4R/run_R04_getDiffMeth_PQLseq_G1.stderr

module load R/4.2.2

cd /data/SBCS-EizaguirreLab/Alice/StickParaBroOff/StickParaOffsBroject/code/R

Rscript R04_getDiffMeth_PQLseq_runInCLUSTER.R
