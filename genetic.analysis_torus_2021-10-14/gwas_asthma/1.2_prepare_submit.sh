#!/bin/bash/


cd $PWD
 

##for oneMCl in Bcell Monocyte NKcell Tcell; do

##for oneMCl in comb comb2; do
oneMCl=comb2
cat motif_split.txt |
while read motifFile; do
   sbatch -q express -p erprp --mem=20G --time=1-00:00:00 -N 1-1 -n 1 --output=slurm_${oneMCl}_${motifFile}.output --wrap "
   module load R;
   R CMD BATCH '--args ${motifFile} ${oneMCl}' 1.2_prepare_torus.R torm_${oneMCl}_${motifFile}.Rout"
   echo ${oneMCl} ${motifFile}
   sleep 1;
done 

##done
