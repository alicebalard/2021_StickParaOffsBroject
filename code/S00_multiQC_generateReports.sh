#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=1G
#$ -l h_rt=240:0:0
#$ -cwd
#$ -j y

# run from conda environment "myenv" where multiqc is loaded
module load anaconda3
conda activate myenv

# Run multiqc on raw reads fastqc quality checks
## multiqc /data/SBCS-EizaguirreLab/Alice/StickParaBroOff/01Quality_Control_fastqc/QCrawReads/

multiqc /data/SBCS-EizaguirreLab/Alice/StickParaBroOff/00Quality_Control_fastqc/QCtrimmedReads/

mv multiqc_data /data/SBCS-EizaguirreLab/Alice/StickParaBroOff/01Trimmed_Reads/.
mv multiqc_report.html /data/SBCS-EizaguirreLab/Alice/StickParaBroOff/01Trimmed_Reads/.

conda deactivate
