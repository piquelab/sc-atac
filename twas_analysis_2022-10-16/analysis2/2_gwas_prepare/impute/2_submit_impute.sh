#!/bin/bash/

cd $PWD

### loop for trait
cat traits.txt | sed -n '4p' | \
while read trait; do 

### loop for split files
cat ./1_missing_snp.outs/${trait}_files.txt | \
while read split_one; do

sbatch -q express -p erprp --mem=40G --time=3-00:00:00 -N 1-1 -n 1 --job-name=impute_${split_one} --output=slurm_impute_${trait}_${split_one}.output --wrap "
module load R;
R CMD BATCH --no-save --no-restore '--args ${trait} ${split_one}' 2_impute.R impute_${trait}_${split_one}.Rout"

echo ${trait} ${split_one}

sleep 1;
done

done
