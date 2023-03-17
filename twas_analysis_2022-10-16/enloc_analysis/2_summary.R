##
library(tidyverse)
library(data.table)
library(devtools)
library(openxlsx)

##
## devtools::install_github("jokamoto97/INTACT", force=T, lib="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
library(INTACT, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")


###
###
gene_enloc <- read.table("enloc.gene.out", header=T)


###
### gwas file 
fn <- "../Asthma_twas/ALOFT_topPIP_twas.txt.gz"
twas <- fread(fn, header=T, data.table=F)


###gwas file
fn <- "../gwas_asthma/asthma_grch37/asthma_summ_ALOFT.txt"
gwas <- fread(fn, header=T, data.table=F)
zscore <- gwas$zscore
names(zscore) <- gwas$id_b37

twas$zscore <- zscore[as.character(twas$genetic_variant)]


### obtain protein coding genes
res2 <- read.xlsx("../4.3_example.outs/ALOFT_topPIP_twas.xlsx")
geneShared <- unique(res2$gene)


###
### Integrate colocalization and TWAS z-score
gene_enloc2 <- gene_enloc%>%filter(Gene%in%geneShared)
twas2 <- twas%>%filter(gene%in%geneShared)%>%dplyr::select(gene, zscore)

DF <- gene_enloc2%>%inner_join(twas2, by=c("Gene"="gene"))
 
res_intact <- intact(GLCP=DF$GLCP, z_vec=DF$zscore)
DF$PCG <- res_intact

###
DF2 <- DF%>%mutate(LFDR=1-PCG)%>%arrange(LFDR)
x <- DF2$LFDR
FDR <- cumsum(x)/1:length(x)
DF2$FDR <- FDR

###
opfn <- "ALOFT_intact.txt"
write.table(DF2, opfn, row.names=F, quote=F)






















