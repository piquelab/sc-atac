#!/bin/bash

module load bedtools 

# cell='Tcell'
for pct0 in 0.05 0.02 0.01
do
   echo ${pct0}
   bedtools intersect -a snp2.bed.gz -b ./4_SNPAnnot.outs/pct_${pct0}/Tcell_Active_peak.bed -c|awk '{print $1, $2, $3, $4}' OFS='\t'|bgzip > ./4_SNPAnnot.outs/pct_${pct0}/Tcell_Active_peak_snp.bed.gz
   ###
done
