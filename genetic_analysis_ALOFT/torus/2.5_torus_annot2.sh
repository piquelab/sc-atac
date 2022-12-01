#!/bin/bash
#SBATCH -q primary
##SBATCH --partition=erprp
#SBATCH --mem=50G
#SBATCH --time=3-00:00:00
#SBATCH -N 1-1
#SBATCH -n 1

module load misc

###########################################################################################################
### comprehensive annotation                                                                            ###
### two kinds of multiple column annotation                                                             ###
### one simplified annotation, union of peaks, union of cell-type motifs, and union of response motifs  ###
###########################################################################################################

####
# outdir="./torus_output/torus_peak_pct_${pct0}_union"

outdir="./torus_output/torus_peak_combineNew"


if [ ! -d ${outdir} ]; then
   mkdir -p ${outdir}
fi

##pct=0.02
##cell=Tcell
##condition=ALOFT

### 1, enrichment analysis
torus -d ./torus_input/ALOFT_eQTL.txt.gz \
  -smap ./torus_input/zzz_snp.map.gz \
  -gmap ./torus_input/zzz_gene.map.gz \
  -annot /nfs/rprdata/julong/sc-atac/genetic.analysis_torus_2021-10-14/SNPannotation/4_SNPAnnot.outs/pct_0.02_combineNew/${ii}_torus.annot.gz \
  -est > ${outdir}/ALOFT_${ii}.est \
  -dump_prior ${outdir}/ALOFT_${ii}_dump.prior   

## 2. qtl analysis  
torus -d ./torus_input/ALOFT_eQTL.txt.gz \
  -smap ./torus_input/zzz_snp.map.gz \
  -gmap ./torus_input/zzz_gene.map.gz \
  -annot /nfs/rprdata/julong/sc-atac/genetic.analysis_torus_2021-10-14/SNPannotation/4_SNPAnnot.outs/pct_0.02_combineNew/${ii}_torus.annot.gz \
  -qtl > ${outdir}/ALOFT_${ii}.rst


