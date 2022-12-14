#!/bin/bash/

condition=Whole_Blood
cat geneList_files.txt | \
while read geneFile;
do
echo ${condition} ${geneFile}
sbatch --export=condition=${condition},geneFile=${geneFile} --output=slurm.${condition}_${geneFile}_out 3_run_dap-g.one.sh
sleep 2;
done
## End 
