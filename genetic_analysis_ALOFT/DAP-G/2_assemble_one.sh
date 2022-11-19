#!/bin/bash
#SBATCH -q express
#SBATCH -p erprp
#SBATCH --mem=20G
#SBATCH --time=2-01:00:00
#SBATCH -N 1-1
#SBATCH -n 1



### Tcell annotation


cat ./geneList/${geneFile} | \
while read ENSG;
do
   ##
   if [ -f ${condition}/${ENSG}.expr -a  -f ${condition}/${ENSG}.cis_snp.bed ]; then 
      echo ${ENSG}
      perl assemble_sbams.pl ALOFT.vcf.gz ${condition}/${ENSG}.cis_snp.bed  ${condition}/${ENSG}.expr ./covar/ALOFT_18_GEPCs.txt > ${condition}/${ENSG}.sbams.dat 
   fi
   ##
done
