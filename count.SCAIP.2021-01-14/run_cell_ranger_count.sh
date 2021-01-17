#!/bin/bash
#SBATCH -q express
#SBATCH --partition=erprp
#SBATCH --mem=120G
#SBATCH --time=2-01:00:00
#SBATCH -N 1-1
#SBATCH -n 10
#SBATCH --output=output_${sample}.log

set -v
set -e 

echo $PWD
cd $PWD

module load bcl2fastq
module load cellranger-atac

## Can use this one as an argument. 
##fastqfolder=../../fastq/
fastqfolder=../fastq/

#transcriptome=/nfs/rprdata/refGenome10x/refdata-cellranger-hg19-1.2.0/
##transcriptome=/nfs/rprdata/refGenome10x/refdata-cellranger-GRCh38-3.0.0/

#refgenome=/nfs/rprdata/refGenome10x/refdata-cellranger-atac-mm10-1.0.0/
#refgenome=/nfs/rprdata/refGenome10x/refdata-cellranger-atac-hg19-1.0.1/
#refgenome=/wsu/home/groups/piquelab/data/refGenome10x/refdata-cellranger-atac-GRCh38-1.2.0/
refgenome=/wsu/home/groups/piquelab/data/refGenome10x/refdata-cellranger-atac-b37-1.2.0

##samplefile=$fastqfolder/outs/input_samplesheet.csv
#find ${fastqfolder} -name 'SCAIP*-ATAC*fastq.gz' | sed 's/.*\///;s/_S.*//' | grep -v 'Undet' | sort | uniq > libList.txt

##cat $samplefile | cut -d, -f2 | grep -v Sample | sort | uniq |\
##find ${fastqfolder} -name 'CC7*fastq.gz' | sed 's/.*CC7/CC7/;s/_S.*.fastq.gz//' | sort | uniq |\
##cat libList.txt | grep -v s[1-5] |
#cat libList.txt | \
#while read sample; 
#do 

fastqs=`find ${fastqfolder} -name "${sample}*fastq.gz" | sed 's/\/SCAIP[^\/]*//g' | sort | uniq`
fastqlist=`echo ${fastqs} | tr ' ' ,`

echo $sample 
echo $fastqlist
 
cellranger-atac count \
      --id=$sample \
      --fastqs=$fastqlist \
      --sample=$sample \
      --reference=$refgenome \
      --localcores=10 --localmem=100 --localvmem=110  



#| qsub -q erprq -l nodes=1:ppn=24 -l mem=120g -N $sample
#     sleep 0.5;
# done

##      --force-cells=5000 \
##      --force-cells=6000 \


## --jobmode=erprq
## qsub -I -q erprq -l nodes=1:ppn=28 -l mem=120g
