#!/bin/bash
# SNP calling on RRBS data using BISCUIT
#$ -N biscuit_SNP_call
#$ -o /data/scratch/btx915/05Biscuit_SNPcall/run_biscuit_SNP_call.stdout
#$ -e /data/scratch/btx915/05Biscuit_SNPcall/run_biscuit_SNP_call.stderr
#$ -cwd
#$ -V
#$ -l h_rt=240:00:00
#$ -l h_vmem=36G

# BISCUIT (BISulfite-seq CUI Toolkit) Version: 1.0.0.dev

cd /data/scratch/btx915/05Biscuit_SNPcall

# index genome
# biscuit index /data/SBCS-EizaguirreLab/Alice/StickParaBroOff/Data/02RefGenome/Gynogen_pchrom_assembly_all.fasta

biscuit align /data/SBCS-EizaguirreLab/Alice/StickParaBroOff/Data/02RefGenome/Gynogen_pchrom_assembly_all.fasta /data/SBCS-EizaguirreLab/Alice/StickParaBroOff/Data/01TrimmedReads_cutadapt/G04986-L1_S114_L007_R1_001_trimmed_cutadapt.fastq.gz
