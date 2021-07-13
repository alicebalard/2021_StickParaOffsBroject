
#!/bin/bash
# produces tab delimited files with content-specific methylation information
# Bismark Methylation Extractor
#$ -N bismark_methylation_extractor
#$ -o /data/scratch/btx915/04Bismark_MethylationExtraction/run_bismark_methylation_extractor.stdout
#$ -e /data/scratch/btx915/04Bismark_MethylationExtraction/run_bismark_methylation_extractor.stderr
#$ -cwd
#$ -V
#$ -l h_rt=240:00:00
#$ -l h_vmem=36G
#$ -t 1-144

module load bismark
module load samtools

DATA_DIR=/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/03Bismark_alignment
OUTDIR=/data/scratch/btx915/04Bismark_MethylationExtraction

cd $DATA_DIR

ls -1 *trimmed_cutadapt_bismark_bt2.sam.gz > list_of_files.txt

# Select the correct line of list of files at each iteration
INPUT_FILE=$(sed -n "${SGE_TASK_ID}p" list_of_files.txt)

bismark_methylation_extractor -o $OUTDIR -s --bedGraph --counts --report $DATA_DIR/$INPUT_FILE


