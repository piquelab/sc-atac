#!/bin/bash/

if [ ! -d Filtered ]; then
   mkdir Filtered
fi

cat libList.txt | \
while read sample; 
do 
   sbatch --export=sample=${sample} --output=./Filtered/slurm.${sample}.output run_filterBam.sh
   sleep 0.5;
done
