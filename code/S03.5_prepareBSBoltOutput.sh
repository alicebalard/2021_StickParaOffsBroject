#!/bin/bash
# BSBolt clean output
#$ -N BSBoltclean
#$ -l h_rt=240:0:0
#$ -e /data/scratch/btx915/BSBolt/MethylationCallingClean/run_BSBolt_clean.stderr
#$ -o /data/scratch/btx915/BSBolt/MethylationCallingClean/run_BSBolt_clean.stdout
#$ -V
#$ -l h_vmem=10G

DIR=/data/scratch/btx915/BSBolt/

cd $DIR

for INPUT_FILE in `ls -1 $DIR/*fastq.gz.bam.sorted.CGmap.gz`
# Separate by C genomic context
do
	gzip -c -d $INPUT_FILE | (awk -F $'\t' '{ if ($4=="CG")  { printf("%s\n", $0); } }' | gzip -c --best > $DIR/MethylationCallingClean/${INPUT_FILE##*/}.CG.map.gz)
	gzip -c -d $INPUT_FILE | (awk -F $'\t' '{ if ($4=="CHG")  { printf("%s\n", $0); } }' | gzip -c --best > $DIR/MethylationCallingClean/${INPUT_FILE##*/}.CHG.map.gz)
	gzip -c -d $INPUT_FILE | (awk -F $'\t' '{ if ($4=="CHH")  { printf("%s\n", $0); } }' | gzip -c --best > $DIR/MethylationCallingClean/${INPUT_FILE##*/}.CHH.map.gz)
	gzip -c -d $INPUT_FILE | (awk -F $'\t' '{ if ($4!="CG" && $4!="CHG" && $4!="CHH")  { printf("%s\n", $0); } }' | gzip -c --best > $DIR/MethylationCallingClean/${INPUT_FILE##*/}.unmap.gz)
done 

# Data moved to: /data/SBCS-EizaguirreLab/Alice/StickParaBroOff/Data/04BSBolt_methCall/BSBolt/MethylationCalling/Methylation_calling_splitted/

