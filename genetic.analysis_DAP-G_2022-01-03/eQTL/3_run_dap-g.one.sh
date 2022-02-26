#!/bin/bash
#SBATCH -q express
#SBATCH -p erprp
#SBATCH --mem=10G
#SBATCH --time=2-01:00:00
#SBATCH -N 1-1
#SBATCH -n 10
##SBATCH --array=1-15770%50

# condition="Bcell_CTRL"

pct0=0.02
outdir="./dap-g_outs/dap-g_pct_${pct0}/${condition}"
if [ ! -d ${outdir} ]; then
   mkdir -p ${outdir}
fi

# i=${SLURM_ARRAY_TASK_ID}
# i=${i}-1

# # ### read into array
# IFS=$'\n' read -d '' -r -a lines < zzz_geneList.txt
# ENSG=${lines[$i]}


cat zzz_geneList.txt | \
while read ENSG;
do
   ##
   if [ -f ${condition}/${ENSG}.sbams.dat ]; then 
      echo ${ENSG}
      /wsu/home/ha/ha21/ha2164/Bin/dap-g.static -d ${condition}/${ENSG}.sbams.dat \
         -p /nfs/rprdata/julong/sc-atac/genetic.analysis_torus_2021-10-14/eQTL/torus_output/torus_peak_pct_${pct0}_response3/${condition}_dump.prior/${ENSG}.prior \
         -t 10 \
         -o  ${outdir}/${ENSG}.out -l ${outdir}/${ENSG}.log 
   fi
   ##
done
