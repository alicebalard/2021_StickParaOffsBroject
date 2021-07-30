#!/bin/bash
# Samtools sort index for BSBolt methylation call
#$ -N samtools
#$ -l h_vmem=36G
#$ -l h_rt=240:0:0
#$ -o /data/scratch/btx915/BSBolt/run_samtools.stdout
#$ -e /data/scratch/btx915/BSBolt/run_samtools.stderr
#$ -V
#$ -t 1-144
#$ -tc 10

module load samtools

BAMDIR=/data/scratch/btx915/BSBolt/Alignments

cd $BAMDIR

# Create the files to loop over:
ls -1 *trimmed_cutadapt.fastq.gz.bam > list_of_files.txt

# Select the correct line of list of files at each iteration
INPUT_FILE=$(sed -n "${SGE_TASK_ID}p" list_of_files.txt)

# sort then index bam file for methylation calling
samtools sort -o ${INPUT_FILE##*/}.sorted.bam $INPUT_FILE

samtools index ${INPUT_FILE##*/}.sorted.bam
