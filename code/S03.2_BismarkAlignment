#!/bin/bash
#Bisulphite sequence alignment of trimmed data on reference gynogen
#$ -N bismark_alignment
#$ -o /data/scratch/btx915/04Bismark_alignment_afterCutadapt/run_bismark_alignment.stdout
#$ -e /data/scratch/btx915/04Bismark_alignment_afterCutadapt/run_bismark_alignment.stderr
#$ -cwd
#$ -V
#$ -l h_rt=240:00:00
#$ -l h_vmem=32G
#$ -t 1-144

module load bismark/0.22.1
module load bowtie2/2.4.1

READS_DIR=/data/scratch/btx915/01TrimmedReads_cutadapt
PATH2GEN=/data/scratch/btx915/03RefGenome
OUTDIR=/data/scratch/btx915/04Bismark_alignment_afterCutadapt

cd $READS_DIR

ls -1 *trimmed_cutadapt.fastq.gz > list_of_files.txt

cd $OUTDIR

INPUT_FILE=$(sed -n "${SGE_TASK_ID}p" $READS_DIR/list_of_files.txt)

for file in `ls $READS_DIR/$INPUT_FILE`
do bismark --genome $PATH2GEN $file; echo $file; done

## Extract mapping efficiency file:
for file in $OUTDIR/*report.txt; do filename=$(echo "${file##*/}") ;  map=$(grep "Mapping efficiency:" $file) ; echo "$filename $map" ; done > Report_mapping_efficiency.txt