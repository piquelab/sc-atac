#!/bin/bash
#SBATCH -q primary
#SBATCH --mem=8G
#SBATCH --time=2-01:00:00
#SBATCH -N 1-1
#SBATCH -n 1
## SBATCH --array=16-20%5 ##1-841%20 ##841


# i=${SLURM_ARRAY_TASK_ID}
# i=${i}-1

# ### read into array
# IFS=$'\n' read -d '' -r -a lines < datasets.txt
# condition=${lines[$i]}
# echo ${condition}

R CMD BATCH 1_prepare.R 1_prepare.Rout 
