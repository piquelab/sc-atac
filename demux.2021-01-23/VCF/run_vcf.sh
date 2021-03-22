#!/bin/bash
#SBATCH -q express
#SBATCH --partition=erprp
#SBATCH --mem=50G
#SBATCH --time=1-01:00:00
#SBATCH -N 1-1
#SBATCH -n 2

set -v
set -e

module load samtools
module load bcftools

cd /nfs/rprdata/julong/sc-atac/demux.2021-01-23/VCF/
echo $PWD


### (2). merge files
bcftools merge `find ./TMP/ -name 'SCAIP6*.bcf'` --threads 10  -Oz -o SCAIP6.merge.vcf.gz 

### (3). extract SNP acoording to position
#user-defined format
bcftools query SCAIP6.merge.vcf.gz -i 'INFO/DP>100' -f '%CHROM\t%POS\n'|bgzip > SCAIP6.posG100.txt.gz
#subset vcf
bcftools view -T <(zcat SCAIP6.posG100.txt.gz) /wsu/home/groups/piquelab/SCAIP/vcf/SCAIP1-6.vcf.gz --threads 10 -Oz -o SCAIP6.posG100.AL.vcf.gz
# #MAF
# bcftools plugin fill-tags SCAIP6.posG100.AL.vcf.gz -- -t 'AN,AC,AF,MAF'|bcftools view --threads 10 -Oz -o SCAIP6.posG100.AF.vcf.gz
# #filter by MAF
# bcftools view SCAIP6.posG100.AF.vcf.gz -i 'INFO/MAF>0.01' --threads 10 -Oz -o SCAIP6.posG100.MAF001.vcf.gz

### (4), extract Batch individuals
bcftools view SCAIP6.posG100.AL.vcf.gz -S ../../sample6.id --threads 2 -Oz -o sample6.posG100.AL.vcf.gz
## MAF
bcftools plugin fill-tags sample6.posG100.AL.vcf.gz -- -t 'AN,AC,AF,MAF'|bcftools view --threads 2 -Oz -o sample6.posG100.AF.vcf.gz 
## filter by MAF
bcftools view sample6.posG100.AF.vcf.gz -i 'INFO/MAF>0.05' --threads 10 -Oz -o sample6.posG100.MAF005.vcf.gz
#bcftools query sample6.posG100.AF.vcf.gz -i 'INFO/MAF>0.1' -f '%CHROM\t%POS\n'|wc -l


### if required to re-order chromosome by alpha-beta, please run the following script.
### But we don't need to re-order chromosome in ATAC project 
### (4). transfer the format required by demuxlet
## we don't need re-order chromosome by alpha-beta in ATAC project.
#reorder VCF lexicographically by chromosome number
#view chromosome order
# bcftools query SCAIP6.posG100.AF.vcf.gz -f '%CHROM\n'|uniq 
# bcftools index SCAIP6.posG100.AF.vcf.gz --threads 10
# bcftools view -r 1,10,11,12,13,14,15,16,17,18,19,2,20,21,22,3,4,5,6,7,8,9 SCAIP6.posG100.AF.vcf.gz --threads 10 -Oz -o SCAIP6.posG100.reordered.vcf.gz
##
# bcftools view -h SCAIP6.posG100.reordered.vcf.gz --threads 10 > my1.vcf.header
# bcftools reheader -h my1.vcf.header SCAIP6.posG100.reordered.vcf.gz --threads 10 -o SCAIP6.posG100.reheader.vcf.gz

##filter by MAF
#bcftools view SCAIP6.posG100.reheader.vcf.gz -i 'INFO/MAF>0.01' --threads 10 -Oz -o SCAIP6.posG100.reheaderMAF.vcf.gz


#zcat SCAIP1-6.merge.posG40.reheader.vcf.gz|head -n 1
#bcftools query SCAIP-ALL.1-6.reordered.vcf.gz -i 'INFO/MAF>0' -f '%CHROM\t%POS\n')|wc -l &
#cat <(bcftools query SCAIP1-6.filtered2.vcf.gz -f '%CHROM\t%POS\t%INFO/AF\n')|wc -l &
#cat <(bcftools query SCAIP1-6.filtered5.vcf.gz -f '%CHROM\t%POS\n' --threads 5)|wc -l &
#bcftools query SCAIP.ALL-1-6.merge.posG4.reordered.vcf.gz -f '%CHROM\t%POS\t%INFO/AN\t%INFO/AC\t%INFO/AF\n') > zzz.posG4.txt
#bcftools view SCAIP-ALL.vcf.gz -h >zzz.header
#cat <(bcftools query SCAIP-ALL.vcf.gz -f '%CHROM\t%POS\t%INFO/AN\t%INFO/AC\n') > zzz.txt
#bcftools query -l SCAIP-ALL.1-6.vcf.gz |head
#zcat SCAIP-ALL.1-6.posG4.txt.gz|wc -l
#bcftools query SCAIP-ALL.1-6.reordered.vcf.gz -f '%CHROM\t%POS\n'|wc -l &
#bcftools query SCAIP-ALL.1-6.AL.vcf.gz -f '%CHROM\t%POS\n'|wc -l &
#bcftools query SCAIP-ALL.1-6.reheader.vcf.gz -i 'INFO/MAF>0.05' -f '%CHROM\t%POS\n'|wc -l &
#bcftools query SCAIP-ALL.1-6.reheader.vcf.gz -f '%CHROM\t%POS\n'|wc -l &
#bcftools query SCAIP-ALL.1-6.merge.vcf.gz -f '%CHROM\t%POS\t%INFO/DP\n'|bgzip > zzz.infor.txt.gz &
#bcftools query SCAIP-ALL.1-6.reheader.vcf.gz -f '%CHROM\t%POS\t%INFO/MAF\n'|bgzip > zzz.infor.txt.gz &
#zcat SCAIP-ALL.1-6.reheader.vcf.gz|sed '2p'
#samtools idxstats bamfile | cut -f 1 |head

#bcftools query SCAIP-ALL.1-6.posG100.reheaderMAF.vcf.gz -f '%CHROM\t%POS\n'|wc -l &
