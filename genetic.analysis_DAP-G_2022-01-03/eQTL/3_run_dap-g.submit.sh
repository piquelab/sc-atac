#!/bin/bash/

cat datasets.txt |sed -n '16,20p'| \
while read condition;
do
   echo ${condition}
   sbatch --export=condition=${condition} --output=slurm.${condition}.out 3_run_dap-g.one.sh
   sleep 1;
done
