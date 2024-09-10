#!/bin/bash

cd $PWD

### loop for traits
cat traits.txt | sort | sed -n '3p' | \
while read trait; do
   nsnp=`zcat ./gwas_PIP2/${trait}.pip.gz |wc -l` 
   echo ${trait} ${nsnp} 
   sbatch -q primary --mem=20G --time=2-10:00:00 -n 1 -N 1-1 --job-name=enloc_${trait} --output=slurm_enloc_${trait}.out --wrap "
   /wsu/home/ha/ha21/ha2164/Bin/fastenloc.static -eqtl fastenloc.eqtl.dtss.vcf.gz -go ./gwas_PIP2/${trait}.pip.gz -total_variants ${nsnp}  -prefix ./enloc_dtss_output/${trait}"
   ###
done
