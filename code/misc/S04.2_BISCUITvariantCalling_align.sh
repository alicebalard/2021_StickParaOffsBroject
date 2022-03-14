#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=20G
#$ -l h_rt=240:0:0
#$ -cwd
#$ -o run_biscuit_Alice.stdout
#$ -e run_biscuit_Alice.stderr
#$ -N biscuit

module load samtools/1.10

OUTDIR=/data/scratch/btx915/05Biscuit_SNPcall
REFGEN=/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/Data/02RefGenome/Gynogen_pchrom_assembly_all.fa
INPUT=/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/Data/01TrimmedReads_cutadapt/G04873-L1_S1_L001_R1_001_trimmed_cutadapt.fastq.gz

## Align
~/bin/biscuit/biscuit align $REFGEN $INPUT | ~/bin/samblaster/samblaster --ignoreUnmated | samtools sort -o $OUTDIR/my_output.bam -O BAM -

samtools index $OUTDIR/my_output.bam

~/bin/biscuit/biscuit pileup -o my_pileup.vcf $REFGEN my_output.bam

bgzip my_pileup.vcf

tabix -p vcf my_pileup.vcf.gz
