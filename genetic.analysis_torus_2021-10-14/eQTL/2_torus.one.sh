#!/bin/bash
#SBATCH -q primary
##SBATCH --partition=erprp
#SBATCH --mem=20G
#SBATCH --time=3-00:00:00
##SBATCH --array=1-5%5
#SBATCH -N 1-1
#SBATCH -n 1

module load misc

# condition_eQTL='Bcell_CTRL'

### create directory
# if [ ! -d "./torus_output/${condition_eQTL}" ]; then
#     mkdir -p "./torus_output/${condition_eQTL}"
# fi

outdir="./torus_output/torus_motif3/${condition_motif}/"
if [ ! -d ${outdir} ]; then
   mkdir -p ${outdir}
fi


###
#i=${SLURM_ARRAY_TASK_ID}
#i=${i}-1
# i=0

# IFS=$'\n' read -d '' -r -a lines < ../SNPannotation/motifList.txt
# motif=${lines[$i]}


### without any annotation
### enrichment
# torus -d ./torus_input/${condition_eQTL}.eQTL.txt.gz \
#   -est > ${outdir}${condition_eQTL}.est

# ### QTL discovery
# torus -d ./torus_input/${condition_eQTL}.eQTL.txt.gz  \
#    -qtl > ${outdir}${condition_eQTL}.egene.rst


### add map information 
### enrichment analysis
# torus -d ./torus_input/${condition_eQTL}.eQTL.txt.gz \
#   -smap ./torus_input/${condition_eQTL}.snp.map.gz \
#   -gmap ./torus_input/${condition_eQTL}.gene.map.gz \
#   -est > ${outdir}${condition_eQTL}.est
#    # # -annot ../SNPannotation/annot_torus/allsnp_${motif}.bed.gz \
#    # -est > ./torus_output/${condition}/${condition}_${motif}.enrichment.est
#    # -qtl > ./torus_output/${condition_eQTL}/${condition_eQTL}.egene.rst
# ### QTL discovery
# torus -d ./torus_input/${condition_eQTL}.eQTL.txt.gz  \
#    -smap ./torus_input/${condition_eQTL}.snp.map.gz \
#    -gmap ./torus_input/${condition_eQTL}.gene.map.gz \
#    -qtl > ${outdir}${condition_eQTL}.egene.rst


### add annotation, position and motif togther  
### enrichment analysis 
# torus -d ./torus_input/${condition_eQTL}.eQTL.txt.gz \
#   -smap ./torus_input/${condition_eQTL}.snp.map.gz \
#   -gmap ./torus_input/${condition_eQTL}.gene.map.gz \
#   -annot ./torus_input/${condition_eQTL}.annot.gz \
#   -est > ${outdir}${condition_eQTL}.est

# ### QTL discovery
# torus -d ./torus_input/${condition_eQTL}.eQTL.txt.gz  \
#    -smap ./torus_input/${condition_eQTL}.snp.map.gz \
#    -gmap ./torus_input/${condition_eQTL}.gene.map.gz \
#    -annot ./torus_input/${condition_eQTL}.annot.gz \
#    -qtl > ${outdir}${condition_eQTL}.egene.rst



### add annotation, position and motif togther  
### enrichment analysis 
torus -d ./torus_input/${condition_eQTL}.eQTL.txt.gz \
  -smap ./torus_input/${condition_eQTL}.snp.map.gz \
  -gmap ./torus_input/${condition_eQTL}.gene.map.gz \
  -annot ./torus_input/annot_positive/zzz_motif.${condition_motif}.annot.gz \
  -est > ${outdir}${condition_eQTL}.est

### QTL discovery
# torus -d ./torus_input/${condition_eQTL}.eQTL.txt.gz  \
#   -smap ./torus_input/${condition_eQTL}.snp.map.gz \
#   -gmap ./torus_input/${condition_eQTL}.gene.map.gz \
#    -annot ./torus_input/zzz_motif.${condition_motif}.annot.gz \
#    -qtl > ${outdir}${condition_eQTL}.egene.rst

  
