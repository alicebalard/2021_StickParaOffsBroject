#!/bin/bash
# Relatedness matrix
#$ -l h_vmem=5G
#$ -l h_rt=1:0:0
#$ -o /data/SBCS-EizaguirreLab/Alice/StickParaBroOff/StickParaOffsBroject/code/logs/S05.2_relatmat.stdout
#$ -e /data/SBCS-EizaguirreLab/Alice/StickParaBroOff/StickParaOffsBroject/code/logs/S05.2_relatmat.stderr
#$ -V

module load bcftools

INPUTDIR=/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/Data/07RelatednessMatrix

## Create a list of files to merge
ls $INPUTDIR/*vcf > $INPUTDIR/list_vcf.txt

## gzip the vcf files (needed for bcftools)
module load htslib/1.19.1

cat $INPUTDIR/list_vcf.txt | while read file
do
    bgzip -c $file > $file.gz
done

## list of zipped files
ls $INPUTDIR/*vcf.gz > $INPUTDIR/list_vcf.gz.txt

## Sort and index files
cat $INPUTDIR/list_vcf.gz.txt | while read file 
do
bcftools sort $file -Oz -o  "${file/input_file/sorted_input_file}" 
bcftools index -t "${file/input_file/sorted_input_file}" 
done

## Merge files
bcftools merge --file-list $INPUTDIR/list_vcf.gz.txt -o $INPUTDIR/merged144samples.vcf.gz

## Rename samples for simplicity
### list old names
ls /data/SBCS-EizaguirreLab/Alice/StickParaBroOff/Data/04BSBolt_methCall/BSBolt/Alignments/*sorted.bam > $INPUTDIR/listbam.tmp

cat $INPUTDIR/listbam.tmp | while read oldname
do                                                                                                                                                                       
## Extract the number after "_S" and before the next "_"
newname=$(echo "$oldname" | sed -n 's/.*_\(S[0-9]*\)_.*/\1/p')
## Append the sample name and file path to the output file, separated by a tab
printf "%s\t%s\n" "$oldname" "$newname" >> $INPUTDIR/renamesamples.txt
done  

## Rename with bcftools
bcftools reheader -s $INPUTDIR/renamesamples.txt -o $INPUTDIR/Merged144samples.vcf.gz $INPUTDIR/merged144samples.vcf.gz
rm $INPUTFIR/merged144samples.vcf.gz

## Convert vcf files to plink ped
module load plink ## v 2.0

## Keep only biallelic variants
plink2 --export ped --vcf $INPUTDIR/Merged144samples.vcf.gz --out $INPUTDIR/Merged --allow-extra-chr --max-alleles 2

## GEMMA requires three files: *.bed, *.bim and *.fam, all with the same prefix.
## The *.bed file should be in the default SNP-major mode (beginning with three bytes).
## One can use the PLINK software to generate binary ped files from standard ped files using the following command:
plink2 --pedmap $INPUTDIR/Merged --make-bed --out $INPUTDIR/Merged --allow-extra-chr

# The .fam file is a plain text file with six columns:
# Family ID: An identifier for the family.
# Individual ID: An identifier for the individual.
# Paternal ID: The ID of the father (0 if unknown).
# Maternal ID: The ID of the mother (0 if unknown).
# Sex: The sex of the individual (1 = male, 2 = female, 0 = unknown/ambiguous).
# Phenotype: The phenotype value (can be a quantitative trait or case/control status).

awk -F, '{print $1, $9, $3, "NA", $13, $12}' /data/SBCS-EizaguirreLab/Alice/StickParaBroOff/StickParaOffsBroject/data/cleanedRawData132fishG1G2.csv | awk 'BEGIN {OFS = " "} {if ($5 == "M") $5 = "1"; else if ($5 == "F") $5 = "2"; print}' |  awk 'BEGIN {OFS = " "} {$4 = "0"; print}' | awk 'BEGIN {OFS = " "} {$3 = ($3 == "NA" ? 0 : $3); print}' | sed '1d' | sort -k2,2V | awk '{$1=$1; print}' OFS='\t' > Merged.fam

## Install GEMMA
GEMMA="/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/StickParaOffsBroject/code/gemma-0.98.5-linux-static-AMD64"

## Create relatedness matrix:
$GEMMA -bfile $INPUTDIR/Merged -gk 1 -o RelatMat132samples

## gk = 1 : centered relatedness matrix

## relatedness matrix in:
$INPUTDIR/output/RelatMat132samples.cXX.txt 
