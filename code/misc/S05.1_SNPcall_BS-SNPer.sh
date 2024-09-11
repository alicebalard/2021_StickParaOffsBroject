#!/bin/bash
# Call SNPs
#$ -l h_vmem=5G
#$ -l h_rt=1:0:0
#$ -o /data/SBCS-EizaguirreLab/Alice/StickParaBroOff/StickParaOffsBroject/code/logs/S05.1_SNPcalls.stdout
#$ -e /data/SBCS-EizaguirreLab/Alice/StickParaBroOff/StickParaOffsBroject/code/logs/S05.1_SNPcalls.stderr
#$ -V
#$ -t 1-144
#$ -tc 10

REF=/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/Data/06GynoGen_functionalAnnotation/Gynogen_pchrom_assembly_all.fasta
BAMSDIR=/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/Data/04BSBolt_methCall/BSBolt/Alignments
OUTDIR=/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/Data/07RelatednessMatrix

## BS-SNPer installed as detailed here: https://github.com/hellbelly/BS-Snper
BSSNPER=/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/StickParaOffsBroject/code/BS-Snper-master/BS-Snper.pl

# Create the files to loop over:
ls -1 $BAMSDIR/*sorted.bam > $OUTDIR/list_of_files.tmp

# Select the correct line of list of files at each iteration
INPUT_FILE=$(sed -n "${SGE_TASK_ID}p" $OUTDIR/list_of_files.tmp)

## Calculate SNPs
perl $BSSNPER --fa $REF $INPUT_FILE --output $OUTDIR/$(basename "$INPUT_FILE")_snp.candidate.out --methcg $OUTDIR/$(basename "$INPUT_FILE")_meth.cg --methchg $OUTDIR/$(basename "$INPUT_FILE")_meth.chg --methchh $OUTDIR/$(basename "$INPUT_FILE")_meth.chh >$OUTDIR/$(basename "$INPUT_FILE")_SNP_BSSNPer_out.vcf 2>$OUTDIR/$(basename "$INPUT_FILE")_SNP.log
