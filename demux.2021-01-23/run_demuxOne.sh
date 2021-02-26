#!/bin/bash
#SBATCH -q primary
##SBATCH --partition=erprp
#SBATCH --mem=300G
#SBATCH --time=6-01:00:00
#SBATCH -N 1-1
#SBATCH -n 1

set -v
set -e

cd /nfs/rprdata/julong/sc-atac/demux.2021-01-23
echo $PWD

module load demuxlet

demuxFolder=demuxOut
vcfFile=./VCF_Filter/sample6.posG100.MAF005.vcf.gz
      
popscle.2021-01-18 dsc-pileup \
  --sam ../count.SCAIP.2021-01-14/Filtered/${sample}.filter_bam.bam \
  --group-list ../count.SCAIP.2021-01-14/${sample}/outs/filtered_peak_bc_matrix/barcodes.tsv \
  --vcf ${vcfFile} \
  --out ${demuxFolder}/${sample}.d > ${demuxFolder}/${sample}.errout.txt;
popscle.2021-01-18 demuxlet --plp ${demuxFolder}/${sample}.d \
  --group-list ../count.SCAIP.2021-01-14/${sample}/outs/filtered_peak_bc_matrix/barcodes.tsv \
  --vcf ${vcfFile} \
  --out ${demuxFolder}/${sample}.out --field GT --alpha 0.0 --alpha 0.5 --doublet-prior 0.1
