#!/bin/bash


###
### create output directory

outdir=./torus_output
if [ ! -d ${outdir} ]; then
   ##
   mkdir -p ${outdir}
fi


###
### submit jobs

for ii in union combine2; do

echo ${ii}

sbatch -q primary --mem=80G --time=14-5:00:00 -n 1 -N 1-1 --job-name=torus_${ii} --output=slurm_torus_${ii}.output --wrap "
  module load misc;
  torus -d ./torus_input/Whole_Blood.eQTL.txt.gz \
    -smap ./torus_input/zzz_snp.map.gz \
    -gmap ./torus_input/zzz_gene.map.gz \
    -annot ./torus_input/zzz_${ii}_torus.annot.gz \
    -est > ${outdir}/Whole_Blood_${ii}.est \
    -dump_prior ${outdir}/Whole_Blood_${ii}_dump.prior"
done

### End

