#!/bin/bash

cd /nfs/rprdata/julong/sc-atac/demux.2021-01-23/VCF_Filter/

if [ ! -d TMP ]; then
   mkdir TMP
fi

cat ../libList.txt |\
while read sample;
do
echo ${sample}
sbatch --export=sample=${sample} --output=slurm.${sample}.out run_mpileup.sh
sleep 1
done
