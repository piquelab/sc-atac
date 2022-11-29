#!/bin/bash


condition=ALOFT

echo ${condition}
dir=./dap-g_outs/dap-g_combineNew_Union/${condition}


cat geneList_files.txt | \
while read geneFile; do
##
echo ${geneFile} 
sbatch --export=dir=${dir},geneFile=${geneFile} --job-name=parse_${geneFile} --output=slurm_parse_${geneFile}.out 4_parse_one.sh

sleep 1;

done





