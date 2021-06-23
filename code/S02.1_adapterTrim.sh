#!/bin/bash
#Trimming RRBS raw illumina data
#$ -N trimgalore
#$ -o run_trimgalore.stdout
#$ -e run_trimgalore.stderr
#$ -cwd
#$ -V
#$ -l h_rt=240:0:0
#$ -l h_vmem=1G

## trimgalore/0.6.5(default) & cutadapt v2.10 on apocrita
module load trimgalore/0.6.5

OUTDIR=/data/scratch/btx915/01Trimmed_Reads/
cd $OUTDIR

DATA_DIR=/data/scratch/btx915/00Illumina_RawReads/

for file in `ls $DATA_DIR/*.fastq.gz`
do
    newname=`basename $file | sed -e "s/.fastq.gz/_trimgalore_cleaned.fastq.gz/"`
    trim_galore -q 20 --phred33 -e 0.1 --length 20 -o $OUTDIR --fastqc_args "--outdir $OUTDIR/FASTQC" --rrbs --gzip --non_directional $file > ./"$newname"
        echo $file
        echo $newname
done

## trimmed reads with extension trimmed.fq.gz (the renaming failed?)

########################
## Trimgalore arguments:  

# -q: Trim low-quality ends from reads in addition to adapter removal. For RRBS samples, quality trimming will be performed first, and adapter trimming is carried in a second round.

# -phred33: Instructs Cutadapt to use ASCII+33 quality scores as Phred scores (Sanger/Illumina 1.9+ encoding) for quality trimming.

# -e: Maximum allowed error rate (no. of errors divided by the length of the matching region)

# --length: Discard reads that became shorter than X because of either quality or adapter trimming.

# --gzip: Compress the output file with gzip

# --fastqc_args "<ARGS>": Passes extra arguments to FastQC

# ‘--rrbs’: for DNA material that was digested with MspI, identifies sequences that were adapter-trimmed and removes another 2 bp from their 3' end. This is to avoid that the filled-in cytosine position close to the second MspI site in a sequence is used for methylation calls. Sequences which were merely trimmed because of poor quality will not be shortened any further.

# ‘--non_directional’: screen adaptertrimmed sequences for the presence of either CAA or CGA at the start of sequences and  clip off the first 2 bases if found. If CAA or CGA are found at the start, no bases will be trimmed off from the 3’ end even if the sequence had some contaminating adapter sequence removed (in this case the sequence read likely originated from either the CTOT or CTOB strand).
