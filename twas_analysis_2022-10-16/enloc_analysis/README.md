# Colocalization analysis

We use fastenloc software to do colocalization analysis. The three steps are as follows,
- We first prepare eqtl file using `perl summarize_dap2enloc.pl -dir dap_rst_dir -vcf ../ALOFT.vcf.gz |gzip - >fastenloc.eqtl.vcf.gz` in the directory `/nfs/rprdata/julong/sc-atac/genetic_analysis_ALOFT/DAP-G/enloc_analysis/`. 
- Prepare gwas files using `sbatch 4_torus_PIP.sh` in the directory `/nfs/rprdata/julong/sc-atac/Asthma_gwas/torus/` followed by `1_prepare.R' to 
- Run `./fastenloc.static -eqtl fastenloc.eqtl.annotation.vcf.gz -go Asthma_gwas_grch37.pip.gz -total_variants 6433121 &`

We also run `2_summmary.R` R script to get the final information, the file named `ALOFT_intact.txt`.


