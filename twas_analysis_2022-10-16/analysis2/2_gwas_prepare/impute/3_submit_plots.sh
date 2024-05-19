#!/bin/bash/

cd $PWD

cat traits.txt | \
while read trait; do

sbatch -q express -p erprp --mem=40G --time=12:00:00 -N 1-1 -n 1 --job-name=plots_${trait} --output=slurm_plots_${trait}.output --wrap "
module load R;
R CMD BATCH --no-save --no-retore '--args ${trait}' 3_plots.R plots_${trait}.Rout"

echo ${trait}
 
done 
