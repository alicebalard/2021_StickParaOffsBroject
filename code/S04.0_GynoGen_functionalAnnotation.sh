#!/bin/bash
# Functional annotation of gynogen genome
#$ -pe smp 1
#$ -l h_vmem=10G
#$ -l h_rt=24:0:0
#$ -o /data/SBCS-EizaguirreLab/Alice/StickParaBroOff/StickParaOffsBroject/code/logs/S04_run_gynofuncal_2024.stdout
#$ -e /data/SBCS-EizaguirreLab/Alice/StickParaBroOff/StickParaOffsBroject/code/logs/S04_run_gynofuncal_2024.stderr
#$ -V

## Source: http://weatherby.genetics.utah.edu/MAKER/wiki/index.php/MAKER_Tutorial_for_WGS_Assembly_and_Annotation_Winter_School_2018#Example_MAKER_Annotation_Project
module load blast+
module load anaconda3/2020.02
module load maker/2.31.9-mpi
## NB interproscan needs java 11, incompatible with R so unload R and reload java 11
module unload R/4.2.2
module load java/11.0.20-openjdk

DIR=/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/Data/06GynoGen_functionalAnnotation
GFF1=$DIR/Gy_allnoM_rd3.maker_apocrita.noseq_corrected.gff
GYNOGENOME1=$DIR/Gynogen_pchrom_assembly_all.fasta
UNIPROTDB=$DIR/uniprot_sprot.fasta

cd $DIR

## Extract proteins from new genome with AGAT. 

## the GFF output from Maker2 needed to be filtered first. Parsing the gff output from maker2 was extremely slow (took >53 hours), but then it took only 3 minutes on a filtered gff. I needed to remove unnecessary gff entries that maker2 had included (match, expressed_sequence_match, protein_match, contig) (Source: https://www.biostars.org/p/9477640/):

# grep -Pv "\tmatch_part\t" $GFF1 | grep -Pv "\tprotein_match\t" | grep -Pv "\texpressed_sequence_match\t" | grep -Pv "\tmatch\t" | grep -Pv "\tcontig\t" > $GFF1.streamlined_for_AGAT.gff

GFF2=$GFF1.streamlined_for_AGAT.gff

### To extract and translate the coding regions:
# agat_sp_extract_sequences.pl -g $GFF1.streamlined_for_AGAT.gff -f $GYNOGENOME1 -p -o $DIR/Gynogen_pchrom_assembly_all_PROTEINS.faa

CODINGFAA=$DIR/Gynogen_pchrom_assembly_all_PROTEINS.faa

## Add functional annotations and meta-data
### Perform external functional analysis. 

# STEP 1: Run blastp command:
# We assigned functional annotation to genes by searching homology against uniprot database (http://www.uniprot.org) using BLASTp alignments. uniprot_sprot.fasta.gz was retrieved on the 02-Mar-2022 at 18:23

# makeblastdb -in $UNIPROTDB -dbtype prot
# blastp -query $DIR/Gynogen_pchrom_assembly_all_PROTEINS.faa -db $UNIPROTDB -num_threads 8 -evalue 1e-6 -max_hsps 1 -max_target_seqs 1 -outfmt 6 -out $DIR/Gynogen_pchrom_assembly_all_PROTEINS.faa.blast.out.blastp

# STEP 2: Run interproscan (improvement Sept 2024 from Charley)
## Install: https://interproscan-docs.readthedocs.io/en/latest/UserDocs.html

## with 10 nodes and 10G for interproscan run one day
## qsub /data/SBCS-EizaguirreLab/Alice/StickParaBroOff/StickParaOffsBroject/code/S04.1_interproscan.sh
 
###############################################################
## Integrate functional annotations into structural annotations
# INTERPROout=$DIR/Gynogen_pchrom_assembly_all_PROTEINS_noast.faa.interproscan.out.tsv

## Input: gff to be modified, coding fasta, blasp output, interproscan tsv output
## Output: a functionally annotated gff3 and bed12

FASTA=Gynogen_pchrom_assembly_all_PROTEINS_noast.faa ## noast=asterisks removed
##FASTA=testInput.100.faa

GFF3=$GFF2
## GFF3=testInput.100.gff3

BLASTPout=$DIR/Gynogen_pchrom_assembly_all_PROTEINS.faa.blast.out.blastp
## BLASTPout=testInput.100.faa.blast.out.blastp

IPRout=$DIR/Gynogen_pchrom_assembly_all_PROTEINS_noast.faa.interproscan.out.tsv ## noast: we removed asterisk for interproscan
##IPRout=testInput.100.faa.interproscan.out.tsv

## Integrate functional annotations into structural annotations.

## 1. add iprscan

##  ipr_update_gff - Takes InterproScan (iptrscan) output and maps domain IDs and GO terms to the Dbxref and Ontology_term attributes in the GFF3 file.
ipr_update_gff $GFF3 $IPRout > $GFF3.funtemp.gff3

## iprscan2gff3 - Takes InterproScan (iprscan) output and generates GFF3 features representing domains. Interesting tier for GBrowse.
iprscan2gff3 $IPRout $GFF3 >> $GFF3.funtemp.gff3

mv $GFF3.funtemp.gff3 $GFF3.fun.gff3

## 2. add blastp

## maker_functional_gff - Maps putative functions identified from BLASTP against UniProt/SwissProt to the MAKER produced GFF3 files in the Note attribute.
maker_functional_gff $UNIPROTDB $BLASTPout $GFF3.fun.gff3 > $GFF3.funtemp.gff3

mv $GFF3.funtemp.gff3 $GFF3.fun.gff3

## if there are 2 occurences of "Note=" in a sentence, remove the first occurence of ";Note=Protein of unknown function" if present
sed '/Note=.*Note=/ {s/;Note=Protein of unknown function//}' $GFF3.fun.gff3 > $GFF3.funtemp.gff3

mv $GFF3.funtemp.gff3 $GFF3.fun.gff3

## Convert gff to bed12 for further analyses
conda activate transdecoder

perl /data/home/btx915/.conda/envs/transdecoder/opt/transdecoder/util/gff3_file_to_bed.pl $GFF3.fun.gff3 > $GFF3.fun.bed12

conda deactivate

