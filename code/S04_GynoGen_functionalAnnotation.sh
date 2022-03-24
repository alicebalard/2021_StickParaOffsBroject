#!/bin/bash
# Functional annotation of gynogen genome
#$ -pe smp 8
#$ -l h_vmem=5G
#$ -l h_rt=240:0:0
#$ -o /data/SBCS-EizaguirreLab/Alice/StickParaBroOff/GIT_StickParaOffsBroject/code/run_gynofuncal.stdout
#$ -e /data/SBCS-EizaguirreLab/Alice/StickParaBroOff/GIT_StickParaOffsBroject/code/run_gynofuncal.stderr
#$ -V

## We assigned functional annotation to genes by searching homology against uniprot database (http://www.uniprot.org) and the reference American genome using BLASTp alignments

## uniprot_sprot.fasta.gz was retrieved on the 02-Mar-2022 at 18:23
## this was blastp -query $GYNOGENOME -db $PROTDB num_threads 8 -evalue 1e-6 -max_hsps 1 -max_target_seqs 1 -outfmt 6 -out $DIR/output.blastp        merged with the reference American genome (Peichel et al. 2017) GCF_016920845.1_GAculeatus_UGA_version5_protein.faa 
## cat GCF_016920845.1_GAculeatus_UGA_version5_protein.faa uniprot_sprot.fasta >> GCF_016920845.1_GAculeatus_UGA_version5_protein_AND_uniprot_sprot_02march2022.faa

# module load anaconda3/2020.02
# conda activate agat
module load blast+

DIR=/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/Data/06GynoGen_functionalAnnotation
GFF1=$DIR/Gy_allnoM_rd3.maker_apocrita.noseq_corrected.gff
GYNOGENOME1=$DIR/Gynogen_pchrom_assembly_all.fasta
PROTDB=$DIR/GCF_016920845.1_GAculeatus_UGA_version5_protein_AND_uniprot_sprot_02march2022.faa
UNIPROTDB=$DIR/uniprot_sprot.fasta

## Extract proteins from new genome with AGAT. 

cd $DIR

## the GFF output from Maker2 needed to be filtered first. Parsing the gff output from maker2 was extremely slow (took >53 hours), but then it took only 3 minutes on a filtered gff. I needed to remove unnecessary gff entries that maker2 had included (match, expressed_sequence_match, protein_match, contig) (Source: https://www.biostars.org/p/9477640/):
# grep -Pv "\tmatch_part\t" $GFF1 | grep -Pv "\tprotein_match\t" | grep -Pv "\texpressed_sequence_match\t" | grep -Pv "\tmatch\t" | grep -Pv "\tcontig\t" > $GFF1.streamlined_for_AGAT.gff

### To extract and translate the coding regions:
# agat_sp_extract_sequences.pl -g $GFF1.streamlined_for_AGAT.gff -f $GYNOGENOME1 -p -o $DIR/Gynogen_pchrom_assembly_all_PROTEINS.faa

## Add functional annotations and meta-data
## Source: http://weatherby.genetics.utah.edu/MAKER/wiki/index.php/MAKER_Tutorial_for_WGS_Assembly_and_Annotation_Winter_School_2018#Add_functional_annotations_and_meta-data

### Perform external functional analysis:
# Run blastp command:

## 1. uniprot + US Peichel stick
# makeblastdb -in $PROTDB -dbtype prot
# blastp -query $DIR/Gynogen_pchrom_assembly_all_PROTEINS.faa -db $PROTDB -num_threads 8 -evalue 1e-6 -max_hsps 1 -max_target_seqs 1 -outfmt 6 -out $DIR/output.blastp

## 2. uniprot only
makeblastdb -in $UNIPROTDB -dbtype prot
blastp -query $DIR/Gynogen_pchrom_assembly_all_PROTEINS.faa -db $UNIPROTDB -num_threads 8 -evalue 1e-6 -max_hsps 1 -max_target_seqs 1 -outfmt 6 -out $DIR/output_onlyuniprot.blastp

## PARSE FUNCTIONAL ANNOTATIONS:
## Chema: I think you can use AGAT to modify your original GFF file incorporating the sequence similarity info into the comments section for each gene (I think you can use this tool from the AGAT package: https://agat.readthedocs.io/en/latest/tools/agat_sq_add_attributes_from_tsv.html).

## agat_sp_manage_functional_annotation.pl:
## The blast against Protein Database (outfmt 6) allows to fill the field/attribute NAME for gene and PRODUCT for mRNA

## 1. uniprot + US Peichel stick
#  TO DO : add PE=to make annotation (like uniprot) https://github.com/NBISweden/AGAT/issues/157

## 2. uniprot only
agat_sp_manage_functional_annotation.pl -f $GFF1.streamlined_for_AGAT.gff -b $DIR/output_onlyuniprot.blastp --db $UNIPROTDB --output $DIR/functiano_uniprot_dir

###########################################################
############################# ERROR works only with uniprot
# module load maker/2.31.9-mpi
## Integrate functional annotations into structural annotations. 
# maker_functional_gff $PROTDB $DIR/output.blastp $NONANNOTGFF > $DIR/Gy_allnoM_rd3.maker_apocrita.noseq_corrected_CURATED.gff
# maker_functional_fasta $PROTDB $DIR/output.blastp $DIR/Gynogen_pchrom_assembly_all_PROTEINS.faa > $DIR/Gynogen_pchrom_assembly_all_PROTEINS_CURATED.faa
