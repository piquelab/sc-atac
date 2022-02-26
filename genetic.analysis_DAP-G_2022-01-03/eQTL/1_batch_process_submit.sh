#!/bin/bash/


cat datasets.txt | \
while read condition;
do
   echo ${condition} 
   sbatch -q express -N 1-1 -n 1 --mem=2G --time=1:00:00 --job-name=${condition} \
      --wrap "perl batch_process.pl -e ./expressions/${condition}.bed.gz -g geno.vcf.gz -c ./covar/${condition}_3.GEPCs.txt -t ${condition}"
done
