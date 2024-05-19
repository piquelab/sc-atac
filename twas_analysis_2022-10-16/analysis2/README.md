# New script for Identification of Asthma risk genes with the combination of TWAS and colocaliation

The working directory is is in the path `/nfs/rprdata/julong/sc-atac/twas_analysis_2022-10-16/analysis2/`.

## gwas data
- `gwas_data` contains the raw gwas data
- `impute` impute the variants in eqtl but not in the gwas data
- `torus` fine-mapping gwas using torus, which will be used in colocalization analysis

 
## colocalization analaysis 
We do colocalization analysis using two cohorts, ALOFT in the folder `3_enloc_aloft` and GTEx in the folder `3_enloc_gtex` 

## TWAS analysis in the folder `4_TWAS_smr` 
We performed the TWAS analysis using SMR approach by integrating the two eQTL cohort (ALOFT cohort and GTEx)  and GBMI_full asthma population.
- `1_ALOFT.R`
- `2_GTEx_wbl.R`
- `3_twas_summary.R`
- `4_INTACT_summ.R`
- Plots and summary table for manuscirpt
 