#!/bin/bash
#SBATCH -q primary
#SBATCH --mem=120G
#SBATCH --time=3-10:00:00
#SBATCH -N 1-1
#SBATCH -n 1

cd $PWD
###submit jobs
/wsu/home/ha/ha21/ha2164/Bin/annovar/annotate_variation.pl \
-filter -dbtype avsnp150 -buildver hg19 SCAIP_NOrs_annovar.input \
/wsu/home/ha/ha21/ha2164/Bin/annovar/humandb/
