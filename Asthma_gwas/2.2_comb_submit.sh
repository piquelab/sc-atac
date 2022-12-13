#!/bin/bash/

cd $PWD

##for ii in comb comb2; do
ii=comb2
cat motif_split.txt | \
while read oneFile; do
  echo ${oneFile}
  sbatch --export=motifFile=${oneFile},comb=${ii} --job-name=torus_${ii}_${oneFile} --output=slurm_${ii}_${oneFile}_torus.out 2.0_comb_one.sh  
  sleep 0.5;
done

##done
