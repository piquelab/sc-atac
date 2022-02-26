#!/bin/bash


cat datasets.txt | \
while read condition;
do
   echo ${condition}
   sbatch -q express -N 1-1 -n 1 --mem=10G --time=12:00:00 --job-name=${condition} \
      --wrap "module load R; R CMD BATCH '--args ${condition}' 1_prepare.R ${condition}_prepare.Rout"
   sleep 0.5
done
