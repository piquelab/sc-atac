#!/bin/bash/

cd $PWD

cat traits.txt | \
while read trait; do
   ###
   sbatch -q primary --mem=40G --time=1-01:00:00 -N 1-1 -n 1 --job-name=torus_${trait} --output=slurm_torus_${trait}.output --wrap "
   module load R;
   R CMD BATCH --no-save --no-retore '--args ${trait}' 1_torus_format.R zzz_torus_${trait}.Rout"
   echo ${trait}
   sleep 1;
done 

  
