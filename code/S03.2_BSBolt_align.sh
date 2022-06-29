#!/bin/bash
# BSBolt align
#$ -N BSBolt
#$ -l h_vmem=36G
#$ -l h_rt=240:0:0
#$ -o /data/scratch/btx915/BSBolt/run_BSBolt_align.stdout
#$ -e /data/scratch/btx915/BSBolt/run_BSBolt_align.stderr
#$ -V
#$ -t 1-144
#$ -tc 10

# BSBolt is activated in a virtual environment
module load bwa/0.7.17
module load python
# module load gcc/10.2.0
source ~/bin/myBSBolt/bin/activate

READS_DIR=/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/Data/01TrimmedReads_cutadapt
OUTDIR=/data/scratch/btx915/BSBolt/Alignments
DATABASE=/data/scratch/btx915/BSBolt/Gynogen_pchrom_assembly_all_DB

cd $OUTDIR

# Step 2: align (https://bsbolt.readthedocs.io/en/latest/align.html)
# Single End Alignment Using Default Commands

# Create the files to loop over:
ls -1 $READS_DIR/*trimmed_cutadapt.fastq.gz > list_of_files.txt

# Select the correct line of list of files at each iteration
INPUT_FILE=$(sed -n "${SGE_TASK_ID}p" list_of_files.txt)

python3 -m bsbolt Align -DB $DATABASE -F1 $INPUT_FILE -O $OUTDIR/BSBoltAlignments_${INPUT_FILE##*/}

deactivate
