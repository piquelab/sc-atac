#!/bin/bash
#SBATCH -q express
#SBATCH -p erprp
#SBATCH --mem=20G
#SBATCH --time=2-01:00:00
#SBATCH -N 1-1
#SBATCH -n 4


###
### different SNP annotation

###
### dtss annotation
outdir="./dap-g_outs/dap-g_combineNew_Union/${condition}"
if [ ! -d ${outdir} ]; then
   mkdir -p ${outdir}
fi

prefix_prior=../torus/torus_output/torus_peak_combineNew/${condition}_Union_dump.prior


cat ./geneList/${geneFile} | \
while read ENSG;
do
   ##
   if [ -f ${condition}/${ENSG}.sbams.dat -a  -f ${prefix_prior}/${ENSG}.prior ]; then 
      echo ${ENSG}
      /wsu/home/ha/ha21/ha2164/Bin/dap-g.static -d ${condition}/${ENSG}.sbams.dat \
         -p ${prefix_prior}/${ENSG}.prior \
         -t 4 \
         --output_all  >  ${outdir}/${ENSG}.out \
         -l ${outdir}/${ENSG}.log 
   fi
   ##
done
