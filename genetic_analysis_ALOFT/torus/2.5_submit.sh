#!/bin/bash

cd $PWD


for ii in  combine2 Union;
do 
 
echo ${ii}
sbatch --export=ii=${ii} --output=slurm_${ii}.out --job-name=torus_${ii} 2.5_torus_annot2.sh

sleep 1;

done
