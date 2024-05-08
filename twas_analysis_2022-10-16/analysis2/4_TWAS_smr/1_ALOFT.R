##
library(tidyverse)
library(data.table)

options(scipen=16)

rm(list=ls())


outdir <- "./twas_smr.outs/"
if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)


### bed 
fn <- "SCAIP_final_bed.gz"
bed <- fread(fn, header=T, data.table=F)
 
chr_pos <- bed$chr_pos_grch38
names(chr_pos) <- bed$genetic_variant



###############################
### using eqtl results 
###############################

###
### FastQTL add grch38 position  
fn <- "/nfs/rprdata/julong/sc-atac/twas_analysis_2022-10-16/ALOFT_results/PC1-18.nominals.eQTL.txt.gz"
fast <- fread(fn, header=F, data.table=F)
names(fast) <- c("gene", "genetic_variant", "DTSS", "pval", "beta")

fast <- fast%>%filter(genetic_variant%in%bed$genetic_variant)%>%
    mutate(chr_pos_grch38=chr_pos[genetic_variant])
fast <- fast%>%mutate(chr=gsub("_.*", "", chr_pos_grch38))%>%filter(chr%in%as.character(1:22))


###
### 
traits <- sort(read.table("traits.txt")$V1)[3:4]
for ( ii in traits){    

### gwas summary data
    
fn <- paste("./gwas_impute/", ii, "_impute.txt.gz", sep="") 
summ <- fread(fn, header=T, data.table=F)
## summ2 <- summ%>%filter(id_b38_0%in%fast$chr_pos_grch38)

## add gwas pval    
pval <- summ$pval
names(pval) <- as.character(summ$id_b38_0)      
fast$pval_gwas <- pval[as.character(fast$chr_pos_grch38)]

## add gwas zscore
zscore <- summ$zscore
names(zscore) <- as.character(summ$id_b38_0)
fast$zscore_gwas <- zscore[as.character(fast$chr_pos_grch38)]
    
### select by eqtl pvalue 
fast2 <- fast%>%drop_na(pval_gwas, zscore_gwas)%>%
    group_by(gene)%>%
    slice_min(order_by=pval, n=1)%>%ungroup()%>%as.data.frame()

### select by pval_gwas    
fast3 <- fast2%>%group_by(gene)%>%
    slice_min(order_by=pval_gwas, n=1, with_ties=F)%>%ungroup()%>%as.data.frame()    
    
###
gfn <- gzfile(paste(outdir, ii, "_aloft_minP_twas.txt.gz", sep=""))
write.table(fast3, gfn, row.names=F, col.names=T, quote=F)

cat(ii, nrow(fast3), "\n")

}    
    
    
    
    

#####################################
### using top PIP  with union   
#######################################


### bed 
## fn <- "SCAIP_final_bed.gz"
## bed <- fread(fn, header=T, data.table=F)

## chr_pos <- bed$chr_pos_grch38
## names(chr_pos) <- bed$genetic_variant



### res <- read.table("./ALOFT_results/dap_topPIP.txt", header=F) ## 19067 to 18185 genes
fn <- "/nfs/rprdata/julong/sc-atac/twas_analysis_2022-10-16/ALOFT_results/aloft_allSNPs_union.txt.gz"
res <- fread(fn, header=F, data.table=F)
names(res) <- c("gene", "genetic_variant", "PIP")

res <- res%>%filter(genetic_variant%in%bed$genetic_variant)%>%
    mutate(chr_pos_grch38=chr_pos[genetic_variant])

res <- res%>%mutate(chr=gsub("_.*", "", chr_pos_grch38))%>%filter(chr%in%as.character(1:22)) 

## res2 <- res%>%group_by(gene)%>%slice_max(order_by=PIP, n=1)%>%ungroup()%>%as.data.frame() ### 17641
 
## res2 <- res2%>%filter(genetic_variant%in%bed$genetic_variant)%>%   ### which drop 883 genes
##     mutate(chr_pos_grch38=chr_pos[genetic_variant])

### top SNP don't find the information in hg38 so dropped
### previously we first assign pval-gwas for all snps for each gene and drop missing value and then select the top PIP
### Probably we need correct 
## res2 <- res2%>%mutate(chr=gsub("_.*", "", chr_pos_grch38))%>%filter(chr%in%as.character(1:22)) ### 16844 gene
## drop <- res2%>%filter(!genetic_variant%in%bed$genetic_variant)
 

###
### 
traits <- sort(read.table("traits.txt")$V1)[3:4]
for ( ii in traits){    
### gwas summary data
fn <- paste("./gwas_impute/", ii, "_impute.txt.gz", sep="") 
summ <- fread(fn, header=T, data.table=F)
## summ2 <- summ%>%filter(id_b38_0%in%res2$chr_pos_grch38)

### Add gwas pval    
pval <- summ$pval
names(pval) <- as.character(summ$id_b38_0)
res$pval_gwas <- pval[as.character(res$chr_pos_grch38)]

### Add gwas zscore    
zscore <- summ$zscore
names(zscore) <- as.character(summ$id_b38_0)
res$zscore_gwas <- zscore[as.character(res$chr_pos_grch38)]

res <- res%>%drop_na(pval_gwas, zscore_gwas)
    
### select by max PIP    
res2 <- res%>%group_by(gene)%>%slice_max(order_by=PIP, n=1)%>%ungroup()%>%as.data.frame() ### 17641

### select by min pval_gwas
res3 <- res2%>%drop_na(pval_gwas)%>%group_by(gene)%>%
    slice_min(order_by=pval_gwas, n=1, with_ties=F)%>%ungroup()%>%as.data.frame()    
    
###
gfn <- gzfile(paste(outdir, ii, "_aloft_topPIP_union_twas.txt.gz", sep=""))
write.table(res3, gfn, row.names=F, col.names=T, quote=F)

cat(ii, nrow(res3), "\n")

}    



#####################################
### top PIP without annotation 
######################################


### bed 
## fn <- "SCAIP_final_bed.gz"
## bed <- fread(fn, header=T, data.table=F)

## chr_pos <- bed$chr_pos_grch38
## names(chr_pos) <- bed$genetic_variant

 
 
### res <- read.table("./ALOFT_results/dap_topPIP.txt", header=F) ## 19067 to 18185 genes
fn <- "/nfs/rprdata/julong/sc-atac/twas_analysis_2022-10-16/ALOFT_results/aloft_allSNPs_dtss.txt.gz"
res <- fread(fn, header=F, data.table=F)
names(res) <- c("gene", "genetic_variant", "PIP")

res <- res%>%filter(genetic_variant%in%bed$genetic_variant)%>%
    mutate(chr_pos_grch38=chr_pos[genetic_variant])

res <- res%>%mutate(chr=gsub("_.*", "", chr_pos_grch38))%>%filter(chr%in%as.character(1:22))

## res2 <- res%>%group_by(gene)%>%slice_max(order_by=PIP, n=1)%>%ungroup()%>%as.data.frame()

## res2 <- res2%>%filter(genetic_variant%in%bed$genetic_variant)%>%
##     mutate(chr_pos_grch38=chr_pos[genetic_variant])
## res2 <- res2%>%mutate(chr=gsub("_.*", "", chr_pos_grch38))%>%filter(chr%in%as.character(1:22))




###
###  
traits <- sort(read.table("traits.txt")$V1)[3:4]
for ( ii in traits){    
### gwas summary data
fn <- paste("./gwas_impute/", ii, "_impute.txt.gz", sep="") 
summ <- fread(fn, header=T, data.table=F)

## summ2 <- summ%>%filter(id_b38_0%in%res2$chr_pos_grch38)

### add pavlue    
pval <- summ$pval
names(pval) <- as.character(summ$id_b38_0)
res$pval_gwas <- pval[as.character(res$chr_pos_grch38)]

### add zscore
zscore <- summ$zscore
names(zscore) <- as.character(summ$id_b38_0)
res$zscore_gwas <- zscore[as.character(res$chr_pos_grch38)]    
    
res <- res%>%drop_na(pval_gwas, zscore_gwas)

### select by max PIP    
res2 <- res%>%group_by(gene)%>%slice_max(order_by=PIP, n=1)%>%ungroup()%>%as.data.frame() ### 17641
### select by min pval_gwas
res3 <- res2%>%drop_na(pval_gwas)%>%group_by(gene)%>%
    slice_min(order_by=pval_gwas, n=1, with_ties=F)%>%ungroup()%>%as.data.frame()        

###
gfn <- gzfile(paste(outdir, ii, "_aloft_topPIP_dtss_twas.txt.gz", sep=""))
write.table(res3, gfn, row.names=F, col.names=T, quote=F)

cat(ii, nrow(res3), "\n")

}    
