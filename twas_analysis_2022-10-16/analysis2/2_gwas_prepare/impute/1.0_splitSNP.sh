#!/bin/bash/

cd $PWD

cat traits.txt | sed -n '4p' | \
while read trait; do
###
  dir=./1_missing_snp.outs/${trait}
  if [ ! -d ${dir} ]; then 
      mkdir -p ${dir}
  fi
  
  split -l 20000 ./1_missing_snp.outs/${trait}_missing.txt -d -a 3 ${dir}/splitGene 
  ls ${dir} | grep split > ./1_missing_snp.outs/${trait}_files.txt
  echo ${trait}
done 
