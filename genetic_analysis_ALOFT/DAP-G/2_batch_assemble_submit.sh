#!/bin/bash/

condition=ALOFT

cat geneList_files.txt | \
while read geneFile;
do
echo ${geneFile}
sbatch --export=condition=${condition},geneFile=${geneFile} --output=slurm_${geneFile}.out --job-name=${condition}_${geneFile} 2_assemble_one.sh
sleep 1;
done


# chmod a+x ${condition}.assemble.cmd

