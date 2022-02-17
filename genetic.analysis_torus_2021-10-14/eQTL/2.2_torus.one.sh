#!/bin/bash
#SBATCH -q primary
##SBATCH --partition=erprp
#SBATCH --mem=20G
#SBATCH --time=3-00:00:00
##SBATCH --array=1-5%5
#SBATCH -N 1-1
#SBATCH -n 1

module load misc

# condition_eQTL='Tcell_PHA-EtOH'
# condition_motif='Bcell_LPS'

### create directory
# if [ ! -d "./torus_output/${condition_eQTL}" ]; then
#     mkdir -p "./torus_output/${condition_eQTL}"
# fi

# outdir="./torus_output/torus_motif2/${condition_motif}/"

###
#i=${SLURM_ARRAY_TASK_ID}
#i=${i}-1
# i=0

# IFS=$'\n' read -d '' -r -a lines < ../SNPannotation/motifList.txt
# motif=${lines[$i]}


### add annotation, position and motif togther  
### enrichment analysis 
# torus -d ./torus_input/${condition_eQTL}.eQTL.txt.gz \
#   -smap ./torus_input/${condition_eQTL}.snp.map.gz \
#   -gmap ./torus_input/${condition_eQTL}.gene.map.gz \
#   -annot ./torus_input/annot_motif/zzz_motif.${condition_motif}.annot.gz \
#   --no_dtss \
#   -est > ${outdir}${condition_eQTL}.est

### QTL discovery
# torus -d ./torus_input/${condition_eQTL}.eQTL.txt.gz  \
#   -smap ./torus_input/${condition_eQTL}.snp.map.gz \
#   -gmap ./torus_input/${condition_eQTL}.gene.map.gz \
#    -annot ./torus_input/zzz_motif.${condition_motif}.annot.gz \
#    -qtl > ${outdir}${condition_eQTL}.egene.rst



####
outdir="./torus_output/torus_peak_pct_${pct0}_response3"
if [ ! -d ${outdir} ]; then
   mkdir -p ${outdir}
fi
cell=${condition_eQTL//_*/}

torus -d ./torus_input/${condition_eQTL}.eQTL.txt.gz \
  -smap ./torus_input/${condition_eQTL}.snp.map.gz \
  -gmap ./torus_input/${condition_eQTL}.gene.map.gz \
  -annot ../SNPannotation/4_SNPAnnot.outs/pct_${pct0}/3_union_${cell}_torus.annot.gz \
  -est > ${outdir}/${condition_eQTL}.est \
  -dump_prior ${outdir}/${condition_eQTL}_dump.prior   
  

