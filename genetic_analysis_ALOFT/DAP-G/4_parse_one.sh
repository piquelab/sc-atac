#!/bin/bash
#SBATCH -q express
#SBATCH -p erprp
#SBATCH --mem=8G
#SBATCH --time=2-01:00:00
#SBATCH -N 1-1
#SBATCH -n 1


cat ./geneList/${geneFile} | \
while read ENSG; do
   ##
   if [ -f ${dir}/${ENSG}.out ]; then
   ##
      echo ${ENSG}
      grep '\[' ${dir}/${ENSG}.out > ${dir}/${ENSG}.model.out
      grep '((' ${dir}/${ENSG}.out > ${dir}/${ENSG}.SNP.out
      grep '{'  ${dir}/${ENSG}.out > ${dir}/${ENSG}.cluster.out
   fi
##
done
