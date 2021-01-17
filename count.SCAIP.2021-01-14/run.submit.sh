#!/bin/bash/

# fastqfolder=../fastq/
# find ${fastqfolder} -name 'SCAIP*-ATAC*fastq.gz' | sed 's/.*\///;s/_S.*//' | grep -v 'Undet' | sort | uniq > libList.txt

cat libList.txt | \
while read sample; 
do 
sbatch --export=sample=${sample} run_cell_ranger_count.sh
sleep 0.5;
done
 
