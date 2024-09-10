# New script for Identification of Asthma risk genes with the combination of TWAS and colocaliation

The working directory is is in the path `/nfs/rprdata/julong/sc-atac/twas_analysis_2022-10-16/analysis2/`.

## Prepare gwas file in the `2_gwas_prepare` 
- `gwas_data` contains the raw gwas data
- `impute` impute the variants that eqtl data has while the gwas data didn't have.  
- `torus` fine-mapping gwas using torus, which will be used in colocalization analysis

 
## colocalization analysis in the `3_enloc_aloft`
We perform colocalization analysis using two cohorts, ALOFT in the folder `3_enloc_aloft` and GTEx in the folder `3_enloc_gtex`. Here we focused ALOFT cohort.
- `2_gwas_rs.R`, We need make transformation for the snp id in the gwas to make these SNPs with the same id to eqtl data.   
- `3_batch_enloc.sh`, Colocalization analysis 

## TWAS analysis in the `4_TWAS_smr` 
We performed the TWAS analysis using SMR approach by integrating the two eQTL cohort (ALOFT cohort and GTEx)  and GBMI_full asthma population.
- `1_ALOFT.R`, TWAS using SMR for ALOFT cohort, the topPIP with annotation and minP
- `2_GTEx_wbl.R`, TWAS using SMR for GTEx Whole blood cohort, the topPIP with annotation and minP
- `3_twas_summary.R`, summarize TWAS results, density, histogram plots and qq plots 
- `4_INTACT_summ.R`, run intact to combine TWAS and colocalization. Also combine all information for output file, including, twas, colocalization, pval and PIP for eqtl and annotation
- `5.1_plots_main.R` and `5.2_plots_supp.R` script for main figures, supp figures and tables in the manuscript
- `5.3_get_genes_table.R` script for supp tables from ALOFT eQTL analysis; `5.3_get_genes_table_gtex.R` for supp table from GTEx WBL eQTL analysis
- `5.4_plots_example_pubs_final.R` script for example genes in main figure from ALOFT eQTL
 

 