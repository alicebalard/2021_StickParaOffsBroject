#!/bin/bash
# Bisulfite convert in-silico the reference genome and index it
#$ -o run_bismark_genome_preparation.stdout
#$ -e run_bismark_genome_preparation.stderr
#$ -cwd
#$ -V
#$ -l h_rt=240:00:00
#$ -l h_vmem=32G

module load bismark/0.22.1
module load bowtie2/2.4.1

PATH2REF=/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/03RefGenome/

cd $PATH2REF

bismark_genome_preparation --bowtie2 --verbose $PATH2REF

