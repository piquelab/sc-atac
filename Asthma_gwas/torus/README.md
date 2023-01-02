## TORUS for asthma gwas 

We first run `./gwas_input/1_make_gwas.R` to add loci column in gwas summary file for torus input". After that, we also need run `1_prepare.R` to keep summary data of genetic variants shared by asthma gwas cohort and annotation file.   
We use the same annotation that were used in eQTL mapping. First run `1_prepare.R` to get the two annotation files by replacing SNP id in ALOFT data using SNP id in the gwas summary data. The two annotation files are as follows, one inclduing single column denoting peaks, also hit by cell-type motifs and by respone motifs, across cell-types and treatments; The other contains 21 columns, representing 1 for peak, 4 for 4 cell-types and 16 for 16 treatments.  
After prepraring all the required input files, we just run the script `2_torus_submit.sh` to do enrichment analysis. The summary script `2_summary.R`.



 


  