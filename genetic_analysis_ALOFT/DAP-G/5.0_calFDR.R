##
library(Matrix)
library(tidyverse)
library(data.table)
## library(clusterProfiler)
## library(org.Hs.eg.db)
## library(annotables)
library(ggplot2)
library(cowplot)
library(RColorBrewer)

rm(list=ls())


###
### calculate FDR



### passing arguments
## args=commandArgs(trailingOnly=T)
## if ( length(args)>0){
##    ##
##    condition <- args[1]
## }else{
##    condition <- "Bcell_CTRL"
## }


###
option <- "dap-g_combineNew_Union"
outdir <- paste("./5_summary.outs/", option, "/", sep="")
if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)

condition <- "ALOFT"




###
### Calculate local FDR for each geen
geneList <- read.table("zzz_geneList.txt")$V1
dap <- NULL
for (ens in geneList){    
###
   dapfn <- paste("./dap-g_outs/", option, "/", condition, "/", ens, ".model.out", sep="")
   if ( file.exists(dapfn)&file.size(dapfn)>0){
      ###    
      resTMP <- read.table(dapfn, fill=T, row.names=NULL, header=F)
      res2 <- resTMP%>%filter(V3==0)%>%mutate(gene=ens)%>%dplyr::select(gene, lfdr=V2)
      dap <- rbind(dap, res2)
  }    
}


### 
dap <- dap%>%mutate(lfdr=as.numeric(lfdr))%>%arrange(lfdr)
x <- dap$lfdr
FDR <- cumsum(x)/1:length(x)
dap$FDR <- FDR
### 
opfn2 <- paste(outdir, condition, "_gene_FDR.txt", sep="")
write.table(dap, opfn2, row.names=F, quote=F, sep="\t")


#####################
## compare egenes ###
#####################

dap1 <- read.table("./5_summary.outs/dap-g_dtss/ALOFT_gene_FDR.txt", header=T)%>%
   filter(FDR<0.1)%>%pull(gene)%>%unique()
##
dap2 <- read.table("./5_summary.outs/dap-g_peak_union/ALOFT_gene_FDR.txt", header=T)%>%
   filter(FDR<0.1)%>%pull(gene)%>%unique()

dap3 <- read.table("./5_summary.outs/dap-g_combineNew_Union/ALOFT_gene_FDR.txt", header=T)%>%
   filter(FDR<0.1)%>%pull(gene)%>%unique()

## ###
## torus1 <- read.table("../torus/3_summary.outs/torus_dtss/summary/ALOFT.rst", header=T)%>%
##    filter(fdr<0.1)%>%pull(gene)%>%unique() 
## torus2 <- read.table("../torus/3_summary.outs/torus_peak_union/summary/ALOFT.rst", header=T)%>%
##    filter(fdr<0.1)%>%pull(gene)%>%unique() 


## length(intersect(torus1, torus2))
## length(torus1)
## length(torus2)

## ##
## ##

## length(intersect(torus1, torus2))
## length(intersect(torus1, dap1))
## length(intersect(torus1, dap2))

## length(intersect(torus2, dap1))
## length(intersect(torus2, dap2))
## length(intersect(dap1, dap2))


length(dap1)
length(dap2)
length(dap3)

length(intersect(dap1, dap2))
length(intersect(dap1, dap3))
length(intersect(dap2, dap3))




### End
