#!/bin/bash
# BSBolt index maker
#$ -N BSBolt
#$ -pe smp 1
#$ -l h_vmem=36G
#$ -l h_rt=240:0:0
#$ -o /data/scratch/btx915/BSBolt/run_BSBolt_index.stdout
#$ -e /data/scratch/btx915/BSBolt/run_BSBolt_index.stderr
#$ -V

# BSBolt is activated in a virtual environment
module load bwa/0.7.17
module load python
source ~/bin/myBSBolt/bin/activate

OUTDIR=/data/scratch/btx915/BSBolt

cd $OUTDIR

# Step 1: RRBS Index, MSPI Cut Format, 40bp Lower Fragment Bound, and 400bp Upper Fragment Bound
python3 -m bsbolt Index -G /data/SBCS-EizaguirreLab/Alice/StickParaBroOff/Data/02RefGenome/Gynogen_pchrom_assembly_all.fa -DB $OUTDIR/Gynogen_pchrom_assembly_all_DB -rrbs -rrbs-cut-format C-CGG -rrbs-lower 40 -rrbs-upper 400

deactivate

## NB: put enough RAM otherwise it does not create the .sa!!!! And blocks at the next step!!!
