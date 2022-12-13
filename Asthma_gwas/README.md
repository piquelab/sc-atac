## TORUS for asthma gwas 

We first run `./gwas_input/1_make_gwas.R` to add loci column in gwas summary file for torus input". After that, we also need run `1_prepare.R` to keep summary data of genetic variants shared by asthma gwas cohort and annotation file.  
In the second step, we run `bash 1.2_prepare_submit.sh` for motifs each by each which call `1.2_prepare_torus.R` to generate the torus format of genetic variants annotation. The script `1.3_submit.sh` which calling 1.3_summary_annot.R summarize the number of SNPs hit by motif.  
We ran `2.2_torus_submit.sh` which calling 2.2_torus_one.sh to do enrichment analysis.
 


  