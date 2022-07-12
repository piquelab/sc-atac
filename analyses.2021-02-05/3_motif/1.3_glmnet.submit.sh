#!/bin/bash
#SBATCH -q primary
#SBATCH --mem=20G
#SBATCH --time=1-01:00:00
#SBATCH -N 1-1
#SBATCH -n 10
#SBATCH --array=12,13,17
##SBATCH --array=1-20%10

module load R

i=${SLURM_ARRAY_TASK_ID}
icell=3
cell="NKcell"
## Bcell, Monocyte, NKcell and Tcell

R CMD BATCH "--args ${i} ${icell} ${cell}" 1.3_glmnet.one.R glmnet_${i}.Rout 
