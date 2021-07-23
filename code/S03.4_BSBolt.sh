#!/bin/bash
# BSBolt align and methylation calling
#$ -N BSBolt
#$ -pe smp 1
#$ -l h_vmem=1G
#$ -l h_rt=240:0:0
#$ -o /data/scratch/btx915/BSBolt/run_BSBolt.stdout
#$ -e /data/scratch/btx915/BSBolt/run_BSBolt.stderr
#$ -V
#$ -t 1-144
#$ -tc 10

# BSBolt is activated by virtualenv
module load python
source ~/bsbolt/bin/activate

READS_DIR=/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/Data/01TrimmedReads_cutadapt
OUTDIR=/data/scratch/btx915/BSBolt

cd $OUTDIR

# Step 1: RRBS Index, MSPI Cut Format, 40bp Lower Fragment Bound, and 400bp Upper Fragment Bound
python3 -m bsbolt Index -G /data/SBCS-EizaguirreLab/Alice/StickParaBroOff/Data/02RefGenome/Gynogen_pchrom_assembly_all.fasta -DB $OUTDIR/Gynogen_pchrom_assembly_all_DB -rrbs -rrbs-cut-format C-CGG -rrbs-lower 40 -rrbs-upper 400

# Step 2: align (https://bsbolt.readthedocs.io/en/latest/align.html)
# Single End Alignment Using Default Commands

# Create the files to loop over:
ls -1 $READS_DIR/*trimmed_cutadapt.fastq.gz > list_of_files.txt

# Select the correct line of list of files at each iteration
INPUT_FILE=$(sed -n "${SGE_TASK_ID}p" list_of_files.txt)

python3 -m bsbolt Align -DB $OUTDIR/Gynogen_pchrom_assembly_all_DB -F1 $INPUT_FILE -O $OUTDIR/BSBolt_alignment_Gynogen

deactivate
