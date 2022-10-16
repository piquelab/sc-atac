#!/bin/bash/

cd $PWD

## loop for each cell-type
# for oneMCl in Bcell Monocyte NKcell
# do
oneMCl=Tcell
  echo ${oneMCl}
  outdir=./torus_output/Asthma_${oneMCl}/
  if [ ! -d ${outdir} ]; then
       mkdir -p ${outdir}
  fi

  cat motif_split.txt | \
  while read oneFile; do
  echo ${oneFile}
  sbatch --export=oneMCl=${oneMCl},motifFile=${oneFile} --job-name=torus_${oneMCl}_${oneFile} --output=slurm_${oneMCl}_${oneFile}_torus.out 2.2_torus_one.sh  
  sleep 1;
  done

# done

### End 


 
     # --wrap "module load misc; 
     #         torus --load_zval -d ./torus_input/Asthma_torus_zval.txt.gz \
     #              -annot ./torus_input/${oneMCl}_union_torus_annot.txt.gz \
     #              -est > ./torus_output/Asthma_${oneMCl}.est "
