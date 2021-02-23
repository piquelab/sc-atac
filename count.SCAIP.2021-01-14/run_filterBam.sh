#!/bin/bash
#SBATCH -q express
#SBATCH --partition=erprp
#SBATCH --mem=120G
#SBATCH --time=2-01:00:00
#SBATCH -N 1-1
#SBATCH -n 5

set -v
set -e 

echo $PWD
cd $PWD

module load samtools/1.11
module load bcftools/1.11


bamFile=./${sample}/outs/possorted_bam.bam
barcodeFile=./${sample}/outs/filtered_peak_bc_matrix/barcodes.tsv
outputFile=./Filtered/${sample}.filter_bam.bam


samtools view \
   -@ 5 \
   -D CB:${barcodeFile} \
   -o ${outputFile} \
   ${bamFile}

samtools index -@ 5 -b ${outputFile}
   

