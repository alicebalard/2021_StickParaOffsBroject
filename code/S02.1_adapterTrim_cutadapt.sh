#!/bin/bash
# Trimming RRBS raw illumina data
#$ -N cutadapt
#$ -o /data/scratch/btx915/01TrimmedReads_cutadapt/run_cutadapt_Alice.stdout
#$ -e /data/scratch/btx915/01TrimmedReads_cutadapt/run_cutadapt_Alice.stderr
#$ -cwd
#$ -V
#$ -l h_rt=240:00:00
#$ -l h_vmem=32G

# to load cutadapt v2.10 we need to load trimgalore, but won't use it
module load trimgalore/0.6.5

OUTDIR=/data/scratch/btx915/01TrimmedReads_cutadapt
cd $OUTDIR

DATA_DIR=/data/archive/archive-SBCS-EizaguirreLab/Alice/StickParaBroOff/00Illumina_RawReads
							       
#cutadapt removes adapter sequences from high-throughput sequencing reads.
# Usage:  cutadapt -a ADAPTER [options] [-o output.fastq] input.fastq

for file in `ls $DATA_DIR/*.fastq.gz`
do
    newname=`basename $file | sed -e "s/.fastq.gz/_trimmed_cutadapt.fastq.gz/"`
    cutadapt -e 0.1 -q 20 -m 20 -O 1 -a NNAGATCGGAAGAGCACAC -a AGATCGGAAGAGCACAC -a ATCGGAAGAGCACAC  -o $OUTDIR/"$newname" $file
    echo $file
    echo $newname
done
