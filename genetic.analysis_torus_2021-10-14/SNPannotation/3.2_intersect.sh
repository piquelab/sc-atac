#!/bin/bash

module load bedtools 

# cell='Tcell'
for cell in Bcell Monocyte NKcell Tcell
do
   echo ${cell}
   bedtools intersect -a snp2.bed.gz -b ./3_SNPAnnot.outs/${cell}/${cell}_Active_peak.bed -c|awk '{print $1, $2, $3, $4}' OFS='\t'|bgzip > ./3_SNPAnnot.outs/${cell}/${cell}_Active_peak_snp.bed.gz 
done
