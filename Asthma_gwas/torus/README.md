## TORUS for asthma gwas 

We first run `./gwas_input/1_make_gwas.R` to add loci column in gwas summary file for torus input". After that, we also need run `1_prepare.R` to keep summary data of genetic variants shared by asthma gwas cohort and annotation file.  
In the second step, we run `bash 1.2_prepare_submit.sh` for motifs each by each which call `1.2_prepare_torus.R` to generate the torus format of genetic variants annotation. The script `1.3_submit.sh` which calling 1.3_summary_annot.R summarize the number of SNPs hit by motif.  
We ran `2.2_torus_submit.sh` which calling 2.2_torus_one.sh to do enrichment analysis. The summary script `2.2_summary.R`. 
We also use the same annotation that were used in eQTL mapping. First run `1_prepare.R` to get the two annotation files, one inclduing single column denoting peaks, also hit by cell-type motifs and by respone motifs, across cell-types and treatments; The other contains 21 columns, representing 1 for peak, 4 for 4 cell-types and 16 for 16 treatments. 
Afther prepraring all the required input files, we just run the script `2_torus_submit.sh`. The summary script `2_summary.R`.



 


  