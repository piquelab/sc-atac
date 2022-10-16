#!/bin/bash
#SBATCH -q express
#SBATCH -p erprp
#SBATCH --mem=10G
#SBATCH --time=24:00:00
#SBATCH -N 1-1
#SBATCH -n 1

module load misc

###
### Just test one Cell-type 
####

outdir=./torus_output/Asthma_${oneMCl}


## torus enrich
cat ./Motif_file/${motifFile} | \
while read motif; do

  if [ -f ./torus_input/Motif_${oneMCl}/${motif}_union_torus.annot.gz ]; then
     echo ${motif}
     torus --load_zval -d ./torus_input/Asthma_torus_zval.txt.gz \
           -annot ./torus_input/Motif_${oneMCl}/${motif}_union_torus.annot.gz \
           -est > ${outdir}/${motif}.est
  fi
  ##
done


##
## -dump_prior ./torus_output/Asthma_Tcell_dump.prior   

