#!/bin/bash

cd /nfs/rprdata/julong/sc-atac/demux.2021-01-23

#ls ../count.SCAIP.2021-01-14 |grep '^SCAIP' |sort >libList.txt
# mkdir ./demuxOut/
cat libList.txt | head -1 | \
while read sample;
do
echo ${sample}
sbatch --export=sample=${sample} --output=slurm.${sample}.out run_demuxOne.sh
sleep 1
done
