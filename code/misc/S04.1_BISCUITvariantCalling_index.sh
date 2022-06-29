
#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=1G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -j y
#$ -N biscuit_index

cd /data/SBCS-EizaguirreLab/Alice/StickParaBroOff/Data/02RefGenome/

~/bin/biscuit/biscuit index Gynogen_pchrom_assembly_all.fa
