#!/bin/bash
#SBATCH -q primary
#SBATCH --mem=20G
#SBATCH --time=2-01:00:00
#SBATCH -N 1-1
#SBATCH -n 10
#SBATCH --array=1-633%20 
##1-841%20 ##841

module load ucscbrowser/2019-05-25

i=${SLURM_ARRAY_TASK_ID}
i=${i}-1

### read into array
IFS=$'\n' read -d '' -r -a lines < motifList2020.txt
motif=${lines[$i]}
echo ${motif}
scanPwmVar ./JASPAR2020_core/${motif}.jaspar -j -base=2.0 /wsu/home/groups/piquelab/data/RefGenome/hg19.2bit snp.bed.gz|\
    awk -v OFS='\t' '$5>10||$8>10 {print $1,$11+1,$11,$1":"$11+1":"$12":"$13}' |sed 's/chr//' |bgzip > ./annot_jaspar2020/allsnp_${motif}.bed.gz
