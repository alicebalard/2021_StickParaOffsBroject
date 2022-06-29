#!/bin/bash
# BSBolt methylation call
#$ -N BSBoltMC
#$ -l h_rt=240:0:0
#$ -o /data/scratch/btx915/BSBolt/run_BSBolt_mc.stdout
#$ -e /data/scratch/btx915/BSBolt/run_BSBolt_mc.stderr
#$ -V
#$ -pe smp 2
#$ -l h_vmem=1G

# BSBolt is activated in a virtual environment
module load python
source ~/bin/myBSBolt/bin/activate

BAMDIR=/data/scratch/btx915/BSBolt/Alignments
OUTDIR=/data/scratch/btx915/BSBolt/MethylationCalling
DATABASE=/data/scratch/btx915/BSBolt/Gynogen_pchrom_assembly_all_DB

cd $OUTDIR

# Step 3: Methylation call. Methylation calling is performed by counting the number of bisulfite converted bases relative to the number of reads observed at each cytonsine. Relative to the reference genome methylation status at a cytosine and guanine can only be called using reads mapped to Watson and Crick strands respectively.

# Loop over files because array failed: needs >=2 nodes per sample!
for INPUT_FILE in `ls -1 $BAMDIR/*trimmed_cutadapt.fastq.gz.bam.sorted.bam`
do
    echo $INPUT_FILE
    # Methylation Calling:
    bsbolt CallMethylation -I $INPUT_FILE -O BSBoltMethCall_$(basename "${INPUT_FILE%.*}") -DB $DATABASE -verbose -t 2 > methylation_stats_$(basename "${INPUT_FILE%.*}").txt
done

deactivate
