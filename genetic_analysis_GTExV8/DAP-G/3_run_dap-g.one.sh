#!/bin/bash
#SBATCH -q express
#SBATCH -p erprp
#SBATCH --mem=40G
#SBATCH --time=2-01:00:00
#SBATCH -N 1-1
#SBATCH -n 4


echo ${condition} ${geneFile}

outdir=./dap-g_outs/Whole_Blood_union
if [ ! -d ${outdir} ]; then
   mkdir -p ${outdir}
fi

###prefix_prior=/nfs/rprdata/Anthony/FUNGEI/Genetics_Resub/Torus_JW/torus_output/annot_hg38_liftover/${condition}_dump.prior

cat ./geneList/${geneFile} | \
while read ENSG;
do
   ##
   ## if [ -f ${condition}/sbams/${ENSG}.sbams.dat ] && [ -f ${prefix_prior}/${ENSG}.prior ]; then 
   if [ -f ${condition}/sbams/${ENSG}.sbams.dat ]; then
      echo ${ENSG}
      /wsu/home/ha/ha21/ha2164/Bin/dap-g.static -d ${condition}/sbams/${ENSG}.sbams.dat \
         -p ./Whole_Blood_union_dump.prior/${ENSG}.prior \
         -t 4 \
         -o  ${outdir}/${ENSG}.rst -l ${outdir}/${ENSG}.log 
    fi
##
done
