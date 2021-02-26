#!/bin/bash
#SBATCH -q primary
##SBATCH --partition=erprp
#SBATCH --mem=50G
#SBATCH --time=1-01:00:00
#SBATCH -N 1-1
#SBATCH -n 1

set -v
set -e

module load samtools/1.11
module load bcftools/1.11

cd /nfs/rprdata/julong/sc-atac/demux.2021-01-23/VCF_Filter/
echo $PWD

### (1). call SNPs for 10 Bam file seprately
samtools mpileup \
  -f /nfs/rprscratch/1Kgenomes/phase2_reference_assembly_sequence/hs37d5.fa.bgz \
  -l <(bcftools query /wsu/home/groups/piquelab/SCAIP/vcf/SCAIP1-6.vcf.gz -f '%CHROM\t%POS\n') \
  /nfs/rprdata/julong/sc-atac/count.SCAIP.2021-01-14/Filtered/${sample}.filter_bam.bam -d 1000000 -g -t DP,AD,ADF,ADR > ./TMP/${sample}.merge.pileup.bcf
##index
bcftools index ./TMP/${sample}.merge.pileup.bcf 
