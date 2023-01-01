## How to get eQTL source files
summary data are from `cp /wsu/home/groups/piquelab/gtex_v8/Data/summary_data/Whole_Blood.allpairs.txt.gz .`   
Obtain SNP list by `zcat Whole_Blood.allpairs.txt.gz |grep -v gene |cut -f 2 |sort|uniq > snpList.txt &`   
 Obtain gene list by `zcat Whole_Blood.allpairs.txt.gz` |grep -v gene |awk '{print $1}' |sort |uniq > geneList.txt &`.
 