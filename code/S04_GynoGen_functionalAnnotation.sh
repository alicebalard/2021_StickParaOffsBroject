#!/bin/bash
# Functional annotation of gynogen genome
#$ -pe smp 1
#$ -l h_vmem=36G
#$ -l h_rt=240:0:0
#$ -o /data/SBCS-EizaguirreLab/Alice/StickParaBroOff/GIT_StickParaOffsBroject/code/run_gynofuncal.stdout
#$ -e /data/SBCS-EizaguirreLab/Alice/StickParaBroOff/GIT_StickParaOffsBroject/code/run_gynofuncal.stderr
#$ -V

## https://link.springer.com/article/10.1186/s12864-016-2523-7
## "Functional annotation of genes was done by searching homology against rice protein sequences of SwissProt (http://www.uniprot.org) using BLASTp alignments with an e-value threshold of 1e-10."

module load anaconda3/2020.02
conda activate agat

## Extract longest isoform from gff made by Maker
# agat_sp_keep_longest_isoform.pl -gff /data/SBCS-EizaguirreLab/Alice/StickParaBroOff/GIT_StickParaOffsBroject/gitignore/bigdata/Gy_allnoM_rd3.maker_apocrita.noseq_corrected.gff -o /data/SBCS-EizaguirreLab/Alice/StickParaBroOff/GIT_StickParaOffsBroject/gitignore/bigdata/Gy_allnoM_rd3.maker_apocrita.noseq_corrected_LONGESTiso.gff

NEWGFF=/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/GIT_StickParaOffsBroject/gitignore/bigdata/Gy_allnoM_rd3.maker_apocrita.noseq_corrected_LONGESTiso.gff

## Extract the coding sequences from the genome, based on Maker structural annotation:
GENOME=/data/SBCS-EizaguirreLab/Gynogen_Genome_Assembly/02_PChrom_Assembly/Gynogen_pchrom_assembly_all.fasta

CODINGFAS=/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/Data/06GynoGen_functionalAnnotation/Gynogen_pchrom_assembly_all.fasta_LONGESTiso_allCodingGenes.fasta

# agat_sp_extract_sequences.pl -g $NEWGFF -f $GENOME -o $CODINGFAS

## Searching homology against the gold standard proteome GCF_016920845.1_GAculeatus_UGA_version5_protein.faa
#module load blast+

REFGENOME=/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/Data/06GynoGen_functionalAnnotation/GCF_016920845.1_GAculeatus_UGA_version5_protein.faa

#makeblastdb -in $REFGENOME -dbtype prot

## BLASTX search protein databases using a translated nucleotide query
BLASTOUT=/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/Data/06GynoGen_functionalAnnotation/blastxGynoAllCodingGenes_vs_GCF_016920845.1.faa.outfmt6

#blastx -query $CODINGFAS -db $REFGENOME -outfmt 6 -evalue 1e-50 -out $BLASTOUT

## 3. Functionaly annotate gff with agat_sp_manage_functional_annotation.pl
## DESCRIPTION
## The script take a gff3 file as input and blast and/or interpro output in order to attach functional annotation to corresponding features within the gff file
## The blast against Protein Database (outfmt 6) allows to fill the field/attribute NAME for gene and PRODUCT for mRNA
agat_sp_manage_functional_annotation.pl -f $NEWGFF -b $BLASTOUT -db $REFGENOME -o /data/SBCS-EizaguirreLab/Alice/StickParaBroOff/Data/06GynoGen_functionalAnnotation/Gy_allnoM_rd3.maker_apocrita.noseq_corrected_functAnnotBlastx.gff

## Change name
cd /data/SBCS-EizaguirreLab/Alice/StickParaBroOff/Data/06GynoGen_functionalAnnotation/Gy_allnoM_rd3.maker_apocrita.noseq_corrected_functAnnotBlastx.gff
mv Gy_allnoM_rd3.maker_apocrita.noseq_corrected_LONGESTiso.gff Gy_allnoM_rd3.maker_apocrita.noseq_corrected_afterAGATcuration.gff
