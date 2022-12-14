#!/bin/bash/

cd $PWD

###################################################################
### torus enrichment analysis for combine2 and Union annotation ###
###################################################################

outdir=./torus_output/Combine
if [ ! -d ${outdir} ]; then
   mkdir -p ${outdir}
fi


## loop for each cell-type
for ii in combine2 Union; do 
###
  sbatch -q primary --mem=50G --time=3:00:00 -N 1-1 -n 1 --job-name=torus_${ii} --output=slurm_${ii}_torus.out --wrap "
  module load misc;               
  torus --load_zval -d ./torus_input/Asthma_torus_zval.txt.gz \
     -annot ./torus_input/Combine/${ii}_torus.annot.gz \
     -est > ${outdir}/Asthma_${ii}.est "  
  sleep 1;
done

### End 


 

