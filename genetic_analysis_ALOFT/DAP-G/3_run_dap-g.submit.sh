#!/bin/bash/

condition=ALOFT
cat geneList_files.txt | \
while read geneFile;
do
echo ${condition} ${geneFile}
sbatch --export=condition=${condition},geneFile=${geneFile} --output=slurm_dtss_${geneFile}.out 3_run_dap-g.one.sh
sleep 1;
done
