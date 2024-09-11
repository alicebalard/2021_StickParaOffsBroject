#!/bin/bash
# Functional annotation of gynogen genome
#$ -pe smp 10
#$ -l h_vmem=10G
#$ -l h_rt=240:0:0
#$ -o /data/SBCS-EizaguirreLab/Alice/StickParaBroOff/StickParaOffsBroject/code/logs/S04_run_gynofuncal.stdout
#$ -e /data/SBCS-EizaguirreLab/Alice/StickParaBroOff/StickParaOffsBroject/code/logs/S04_run_gynofuncal.stderr
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

NEWGFF=$DIR/Gy_allnoM_rd3.maker_apocrita.noseq_corrected.gff.streamlined_for_AGAT.gff

### To extract and translate the coding regions:
# agat_sp_extract_sequences.pl -g $GFF1.streamlined_for_AGAT.gff -f $GYNOGENOME1 -p -o $DIR/Gynogen_pchrom_assembly_all_PROTEINS.faa

## Add functional annotations and meta-data
### Perform external functional analysis. 

# STEP 1: Run blastp command:
# We assigned functional annotation to genes by searching homology against uniprot database (http://www.uniprot.org) using BLASTp alignments. uniprot_sprot.fasta.gz was retrieved on the 02-Mar-2022 at 18:23

# makeblastdb -in $UNIPROTDB -dbtype prot
# blastp -query $DIR/Gynogen_pchrom_assembly_all_PROTEINS.faa -db $UNIPROTDB -num_threads 8 -evalue 1e-6 -max_hsps 1 -max_target_seqs 1 -outfmt 6 -out $DIR/output_onlyuniprot.blastp

BLASTPout=$DIR/output_onlyuniprot.blastp

# STEP 2: Run interproscan (improvement Sept 2024 from Charley)
## Install: https://interproscan-docs.readthedocs.io/en/latest/UserDocs.html

## Gene3D (4.3.0) : Structural assignment for whole genes and genomes using the CATH domain structure database.
## PANTHER (15.0) : The PANTHER (Protein ANalysis THrough Evolutionary Relationships) Classification System is a unique resource that classifies genes by their functions, using published scientific experimental evidence and evolutionary relationships to predict function even in the absence of direct experimental evidence.
## Pfam (34.0) : A large collection of protein families, each represented by multiple sequence alignments and hidden Markov models (HMMs)
## SFLD (4) : SFLD is a database of protein families based on hidden Markov models (HMMs).
## SUPERFAMILY (1.75) : SUPERFAMILY is a database of structural and functional annotations for all proteins and genomes.
## TIGRFAM (15.0) : TIGRFAMs are protein families based on hidden Markov models (HMMs).

INTERPROSCAN=./my_interproscan/interproscan-5.56-89.0/interproscan.sh

echo "ncpu:"
echo $NSLOTS

# Set variables
INPUT_FILE="Gynogen_pchrom_assembly_all_PROTEINS_noast.faa"
OUTPUT_PREFIX="Gynogen_pchrom_assembly_all_PROTEINS_noast.faa.iprscan"

echo "Run InterProScan..."
# Run InterProScan
$INTERPROSCAN -appl Gene3D,PANTHER,pfam,PIRSR,SFLD,SUPERFAMILY,TIGRFAM \
  -i $INPUT_FILE \
  -f TSV,GFF3 \
  -t p \ 
  -o $OUTPUT_PREFIX \
  -cpu $NSLOTS \
  --goterms \
  --pathways \
  -dp -iprlookup
echo "InterProScan analysis complete."

IPRout=$DIR/Gynogen_pchrom_assembly_all_PROTEINS_noast.faa.iprscan

## Integrate functional annotations into structural annotations. 
##  ipr_update_gff - Takes InterproScan (iptrscan) output and maps domain IDs and GO terms to the Dbxref and Ontology_term attributes in the GFF3 file.
ipr_update_gff $NEWGFF $IPRout > $DIR/temp.gff

## iprscan2gff3 - Takes InterproScan (iprscan) output and generates GFF3 features representing domains. Interesting tier for GBrowse.
iprscan2gff3 $IPRout $NEWGFF >> $DIR/temp.gff

## maker_functional_gff - Maps putative functions identified from BLASTP against UniProt/SwissProt to the MAKER produced GFF3 files in the Note attribute.
maker_functional_gff $UNIPROTDB $BLASTPout $DIR/temp.gff > Gy_allnoM_rd3.maker_apocrita.noseq_corrected.gff.streamlined_for_AGAT.CURATED.gff

## Remove double "Note"
sed -i 's/;Note=Protein of unknown function;Note/;Note/g' Gy_allnoM_rd3.maker_apocrita.noseq_corrected.gff.streamlined_for_AGAT.CURATED.gff

## Convert gff to bed12 for further analyses
conda activate transdecoder

perl /data/home/btx915/.conda/envs/transdecoder/opt/transdecoder/util/gff3_file_to_bed.pl $DIR/Gy_allnoM_rd3.maker_apocrita.noseq_corrected.gff.streamlined_for_AGAT.CURATED.gff > $DIR/Gy_allnoM_rd3.maker_apocrita.noseq_corrected.gff.streamlined_for_AGAT.CURATED.transdec.bed12

conda deactivate
