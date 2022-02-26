#!/bin/bash
#SBATCH -q express
#SBATCH -p erprp
#SBATCH --mem=2G
#SBATCH --time=2-01:00:00
#SBATCH -N 1-1
#SBATCH -n 1
#SBATCH --array=16-20


i=${SLURM_ARRAY_TASK_ID}
i=${i}-1

# # ### read into array
IFS=$'\n' read -d '' -r -a lines <  datasets.txt

condition=${lines[$i]}


echo ${condition}
dir="./dap-g_outs/dap-g_pct_0.02/${condition}"

##
cat zzz_geneList.txt | \
while read ENSG;
do
##
   if [ -f ${dir}/${ENSG}.out ]; then
   ##
      grep '\[' ${dir}/${ENSG}.out > ${dir}/${ENSG}.model.out
      grep '((' ${dir}/${ENSG}.out > ${dir}/${ENSG}.SNP.out
      grep '{'  ${dir}/${ENSG}.out > ${dir}/${ENSG}.cluster.out
   fi
##
done





