#!/bin/bash
#SBATCH -q express
#SBATCH -p erprp
#SBATCH --mem=10G
#SBATCH --time=24:00:00
#SBATCH -N 1-1
#SBATCH -n 1

module load misc

###
### TORUS enrichment analysis for each motif combine all the cell-types

outdir=./torus_output/${comb}
if [ ! -d ${outdir} ]; then
   ##
   mkdir -p ${outdir}
fi


## torus enrich
cat ./Motif_file/${motifFile} | \
while read motif; do

  if [ -f ./torus_input/Motif_${comb}/${motif}_torus.annot.gz ]; then
     echo ${motif}
     torus --load_zval -d ./torus_input/Asthma_torus_zval.txt.gz \
           -annot ./torus_input/Motif_${comb}/${motif}_torus.annot.gz \
           -est > ${outdir}/${motif}.est
  fi
  ##
done


##
## -dump_prior ./torus_output/Asthma_Tcell_dump.prior   

