#!/bin/bash/


### 
# cat datasets_motif.txt | \
# while read condition_motif;
# do

   ###
   # echo ${condition_motif}
   # cat datasets.txt | \
   # while read condition_eqtl;
   # do
   #    echo ${condition_eqtl}
   #    # sbatch --export=condition_eQTL=${condition_eqtl},condition_motif=${condition_motif} --output=slurm.${condition_motif}.${condition_eqtl}.out 2_torus.one.sh
   #    sbatch --export=condition_eQTL=${condition_eqtl},cell=${cell} --output=slurm.${condition_eqtl}.${cell}.out 2.2_torus.one.sh 
   #    sleep 1;
   # done
   ###
# done



###
###
# for th in 0.05 0.02 0.01
# do
th=0.01
cat datasets.txt | sed -n '16,20p'| \
while read condition_eqtl;
do
  echo ${condition_eqtl} ${th}
  # sbatch --export=condition_xeQTL=${condition_eqtl},condition_motif=${condition_motif} --output=slurm.${condition_motif}.${condition_eqtl}.out 2_torus.one.sh
  sbatch --export=condition_eQTL=${condition_eqtl},pct0=${th} --output=slurm.${th}.${condition_eqtl}.out 2.2_torus.one.sh 
  sleep 1;
done

# done
### enrichment for cell type separately
# for cell in Bcell Monocyte NKcell Tcell;
# do
#    cat datasets.txt | \
#    while read condition_eqtl;
#    do
#       echo ${condition_eqtl}
#       # sbatch --export=condition_eQTL=${condition_eqtl},condition_motif=${condition_motif} --output=slurm.${condition_motif}.${condition_eqtl}.out 2_torus.one.sh
#       sbatch --export=condition_eQTL=${condition_eqtl},cell=${cell} --output=slurm.${condition_eqtl}.${cell}.out 2.2_torus.one.sh 
#       sleep 1;
#    done
#    ###
# done
