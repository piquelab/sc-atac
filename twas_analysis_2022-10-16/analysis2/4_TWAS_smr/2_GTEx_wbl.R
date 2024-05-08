##
library(tidyverse)
library(data.table)

rm(list=ls())


### passing argument
## args=commandArgs(trailingOnly=T)
## if (length(args)>0){
##    condition <- args[1]
##    dir <- args[2]
##  }else{
##    condition <- "Bcell_CTRL"
##    dir <- "SCAIP_results"
## }

 

###
outdir <- "./twas_smr.outs/"
if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)



###
### 1. FastQTL add grch38 position   
fn <- "/nfs/rprdata/julong/sc-atac/twas_analysis_2022-10-16/Whole_Blood_GTEx_v8_results/Whole_Blood.fastqtl.gz"
fast <- fread(fn, header=F, data.table=F)
names(fast) <- c("gene", "genetic_variant", "DTSS", "pval", "beta", "se")

## fast2 <- fast%>%group_by(gene)%>%slice_min(order_by=pval, n=1)%>%ungroup()%>%as.data.frame()
infor <- str_split(gsub("chr", "", fast$genetic_variant), "_", simplify=T)%>%as.data.frame()
infor <- infor[,1:2]
names(infor) <- c("chr", "pos")
infor <- infor%>%mutate(chr_pos_grch38=paste(chr, pos, sep="_"))

fast <- cbind(fast, infor)%>%filter(chr%in%as.character(1:22))


### 
traits <- sort(read.table("traits.txt")$V1)[3:4]
for ( ii in traits){    
### gwas summary data
fn <- paste("./gwas_impute/", ii, "_impute.txt.gz", sep="") 
summ <- fread(fn, header=T, data.table=F)
summ2 <- summ%>%filter(id_b38_0%in%fast$chr_pos_grch38)

###    
### add gwas pval    
pval <- summ2$pval
names(pval) <- as.character(summ2$id_b38_0)      
fast$pval_gwas <- pval[as.character(fast$chr_pos_grch38)]
    
### add zscore
zscore <- summ2$zscore
names(zscore) <- as.character(summ2$id_b38_0)
fast$zscore_gwas <- zscore[as.character(fast$chr_pos_grch38)]

    
###
### select by eqtl pvalue 
fast2 <- fast%>%drop_na(pval_gwas, zscore_gwas)%>%
    group_by(gene)%>%
    slice_min(order_by=pval, n=1)%>%ungroup()%>%as.data.frame()

### select by pval_gwas    
fast3 <- fast2%>%group_by(gene)%>%
    slice_min(order_by=pval_gwas, n=1, with_ties=F)%>%ungroup()%>%as.data.frame()
    
fast3 <- fast3%>%mutate(gene=gsub("\\..*", "", gene))    
    
###
gfn <- gzfile(paste(outdir, ii, "_gtex_minP_twas.txt.gz", sep=""))
write.table(fast3, gfn, row.names=F, col.names=T, quote=F)

cat(ii, nrow(fast3), "\n")

}    




##########################################
### topPIP  with union annotation 
#########################################

###
### 1. PIP results  add grch38 position   
fn <- "/nfs/rprdata/julong/sc-atac/twas_analysis_2022-10-16/Whole_Blood_GTEx_v8_results/Whole_Blood_allSNPs_union.txt.gz"
res <- fread(fn, header=F, data.table=F)
names(res) <- c("gene", "genetic_variant", "PIP")

## res2 <- res%>%group_by(gene)%>%slice_max(order_by=PIP, n=1)%>%ungroup()%>%as.data.frame()
infor <- str_split(gsub("chr", "", res$genetic_variant), "_", simplify=T)%>%as.data.frame()
infor <- infor[,1:2]
names(infor) <- c("chr", "pos")
infor <- infor%>%mutate(chr_pos_grch38=paste(chr, pos, sep="_"))

res <- cbind(res, infor)%>%filter(chr%in%as.character(1:22))



###
### 
traits <- sort(read.table("traits.txt")$V1)[3:4]
for ( ii in traits[2]){    
### gwas summary data
fn <- paste("./gwas_impute/", ii, "_impute.txt.gz", sep="") 
summ <- fread(fn, header=T, data.table=F)
summ2 <- summ%>%filter(id_b38_0%in%res$chr_pos_grch38)

####    
###  add pval    
pval <- summ2$pval
names(pval) <- as.character(summ2$id_b38_0)    
res$pval_gwas <- pval[as.character(res$chr_pos_grch38)]

## add zscore
zscore <- summ2$zscore
names(zscore) <- as.character(summ2$id_b38_0)
res$zscore_gwas <- zscore[as.character(res$chr_pos_grch38)]


###
### select variants by maxmium PIP for each gene
res2 <- res%>%group_by(gene)%>%slice_max(order_by=PIP, n=1)%>%ungroup()%>%as.data.frame()

### select variants by minmum pval_gwas for each gene    
res3 <- res2%>%drop_na(pval_gwas)%>%group_by(gene)%>%
    slice_min(order_by=pval_gwas, n=1, with_ties=F)%>%ungroup()%>%as.data.frame()    

res3 <- res3%>%mutate(gene=gsub("\\..*", "", gene))    
    
###
gfn <- gzfile(paste(outdir, ii, "_gtex_topPIP_union_twas.txt.gz", sep=""))
write.table(res3, gfn, row.names=F, col.names=T, quote=F)

cat(ii, nrow(res3), "\n")

}    




########################################
### topPIP with dtss annotation  
########################################

###
### 1. PIP results  add grch38 position    
fn <- "/nfs/rprdata/julong/sc-atac/twas_analysis_2022-10-16/Whole_Blood_GTEx_v8_results/Whole_Blood_allSNPs_dtss.txt.gz"
res <- fread(fn, header=F, data.table=F)
names(res) <- c("gene", "genetic_variant", "PIP")

## res2 <- res%>%group_by(gene)%>%slice_max(order_by=PIP, n=1)%>%ungroup()%>%as.data.frame()
infor <- str_split(gsub("chr", "", res$genetic_variant), "_", simplify=T)%>%as.data.frame()
infor <- infor[,1:2]
names(infor) <- c("chr", "pos")
infor <- infor%>%mutate(chr_pos_grch38=paste(chr, pos, sep="_"))

res <- cbind(res, infor)%>%filter(chr%in%as.character(1:22))



###
### 
traits <- sort(read.table("traits.txt")$V1)[3:4]
for ( ii in traits){    
### gwas summary data
fn <- paste("./gwas_impute/", ii, "_impute.txt.gz", sep="") 
summ <- fread(fn, header=T, data.table=F)
summ2 <- summ%>%filter(id_b38_0%in%res$chr_pos_grch38)

### add pvalue     
pval <- summ2$pval
names(pval) <- as.character(summ2$id_b38_0)    
res$pval_gwas <- pval[as.character(res$chr_pos_grch38)]

### add zscore
zscore <- summ2$zscore
names(zscore) <- as.character(summ2$id_b38_0)
res$zscore_gwas <- zscore[as.character(res$chr_pos_grch38)]


###
### select variants by maxmium PIP for each gene    
res2 <- res%>%group_by(gene)%>%slice_max(order_by=PIP, n=1)%>%ungroup()%>%as.data.frame()

### select variants by minmum pval_gwas for each gene    
res3 <- res2%>%drop_na(pval_gwas)%>%group_by(gene)%>%
    slice_min(order_by=pval_gwas, n=1, with_ties=F)%>%ungroup()%>%as.data.frame()    

res3 <- res3%>%mutate(gene=gsub("\\..*", "", gene))    
    
###
gfn <- gzfile(paste(outdir, ii, "_gtex_topPIP_dtss_twas.txt.gz", sep=""))
write.table(res3, gfn, row.names=F, col.names=T, quote=F)

cat(ii, nrow(res3), "\n")

}    

 
###
### End
