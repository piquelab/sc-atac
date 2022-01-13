#!/bin/bash/

cat datasets_motif.txt | \
while read condition_motif;
do

   ###
   echo ${condition_motif}
   cat datasets.txt | \
   while read condition_eqtl;
   do
      echo ${condition_eqtl}
      sbatch --export=condition_eQTL=${condition_eqtl},condition_motif=${condition_motif} --output=slurm.${condition_motif}.${condition_eqtl}.out 2_torus.one.sh 
      sleep 1;
   done
   ###

done


###
#    if [ ! -d './torus_output/${condition_eQTL}' ]; then
#       mkdir -p './torus_output/${condition_eQTL}'
#    fi
# ###
#    cat datasets_motif.txt|\
#    while read condition_motif;
#    do
#       echo ${condition_eQTL}  ${condition_motif}
#    done
## sed -n '14p' | \
## '3p;5p;7p;10,12p'| \
