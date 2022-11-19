### script for running DAP-G

### 1st step, generate ENSG.expr and ENSG.cis_snp.bed files
```ruby
perl batch_process.pl -e ./expressions/ALOFT.bed.gz -g ALOFT.vcf.gz -c ./covar/ALOFT_18_GEPCs.txt -t ALOFT
```

### 2rd step, assemble expression, covariates and SNPs into ENSG.sbams.dat file
```ruby
bash 2_batch_assemble_submit.sh
```

### 3rd step
#### 3.1
```ruby
cd geneList
split zzz_geneList.txt -l 1000 -d -a 3 splitGene
```

#### 3.2, submit jobs
```ruby
bash 3_run_dap-g.submit.sh
```
