#!/bin/bash
# Functional annotation of gynogen genome interproscan
#$ -pe smp 10
#$ -l h_vmem=10G
#$ -l h_rt=240:0:0
#$ -o /data/SBCS-EizaguirreLab/Alice/StickParaBroOff/StickParaOffsBroject/code/logs/S04.1interproscan.stdout
#$ -e /data/SBCS-EizaguirreLab/Alice/StickParaBroOff/StickParaOffsBroject/code/logs/S04.1interproscan.stderr
#$ -V

## with 10 nodes and 10G for interproscan run one day

module load blast+
module load anaconda3/2020.02
module load maker/2.31.9-mpi
## NB interproscan needs java 11, incompatible with R so unload R and reload java 11
module unload R/4.2.2
module load java/11.0.20-openjdk

DIR=/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/Data/06GynoGen_functionalAnnotation

cd $DIR

INTERPROSCAN=./my_interproscan/interproscan-5.56-89.0/interproscan.sh

echo "ncpu:"
echo $NSLOTS

# Set variables
INPUT_FILE="Gynogen_pchrom_assembly_all_PROTEINS_noast.faa"

echo "Run InterProScan..."

$INTERPROSCAN -appl Gene3D,PANTHER,pfam,PIRSR,SFLD,SUPERFAMILY,TIGRFAM  -i $INPUT_FILE -f TSV,GFF3,XML -t p -cpu $NSLOTS --goterms --pathways -dp -iprlookup -b $INPUT_FILE.interproscan.out

echo "InterProScan analysis complete."

## Gene3D (4.3.0) : Structural assignment for whole genes and genomes using the CATH domain structure database.
## PANTHER (15.0) : The PANTHER (Protein ANalysis THrough Evolutionary Relationships) Classification System is a unique resource that classifies genes by their functions, using published scientific experimental evidence and evolutionary relationships to predict function even in the absence of direct experimental evidence.
## Pfam (34.0) : A large collection of protein families, each represented by multiple sequence alignments and hidden Markov models (HMMs)
## SFLD (4) : SFLD is a database of protein families based on hidden Markov models (HMMs).
## SUPERFAMILY (1.75) : SUPERFAMILY is a database of structural and functional annotations for all proteins and genomes.
## TIGRFAM (15.0) : TIGRFAMs are protein families based on hidden Markov models (HMMs).
