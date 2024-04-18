# Identify asthma risk genes with combination of colocalization and TWAS



We first downloaded the data in the folder `./gwas_asthma_GBMI/`, including three ancestry population summary data, African, European and all the populations. The folder `./gwas_asthma_BGMI/` include two types of analysis,
- In the subfolder `./gwas_data/`, we find the representative SNPs for each gene from (ALOFT and GTEx) in fastqtl and dap-g that don't appear in three gwas datasets respectively and return chr_pos_grch38 position. For these missing SNPs, we use the closest SNP within 10 kb region as subsitition. In final we generated a new gwas data called imputation files, which will be used for twas-smr. 
- The subfolder `./torus/`, using torus fine-mapping gwas variants, which be the input data of colocalization.     

The folder `analysis2` for new gwas data, including the following procedures,