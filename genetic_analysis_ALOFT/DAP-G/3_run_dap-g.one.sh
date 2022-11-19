#!/bin/bash
#SBATCH -q express
#SBATCH -p erprp
#SBATCH --mem=20G
#SBATCH --time=2-01:00:00
#SBATCH -N 1-1
#SBATCH -n 4



### Tcell annotation
# outdir="./dap-g_outs/dap-g_peak_union/${condition}"
# if [ ! -d ${outdir} ]; then
#    mkdir -p ${outdir}
# fi

# prefix_prior=../torus/torus_output/torus_peak_union/${condition}_dump.prior


###
### dtss annotation
outdir="./dap-g_outs/dap-g_dtss/${condition}"
if [ ! -d ${outdir} ]; then
   mkdir -p ${outdir}
fi

prefix_prior=../torus/torus_output/torus_dtss/${condition}_dump.prior


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
