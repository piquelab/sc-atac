#!/bin/bash/


cd $PWD
 
# oneMCl="Tcell"

# for oneMCl in Bcell Monocyte NKcell Tcell; do
oneMCl="comb"
cat motif_split.txt |
while read motifFile; do
   sbatch -q express -p erprp --mem=8G --time=1-00:00:00 -N 1-1 -n 1 --output=slurm_${oneMCl}_${motifFile}.output --wrap "
   module load R;
   R CMD BATCH '--args ${motifFile} ${oneMCl}' 1.3_summary_annot.R torm_${oneMCl}_${motifFile}.Rout"
   echo ${oneMCl} ${motifFile}
##   sleep 1;
done 

###done
