##
library(tidyverse)
library(data.table)
options(scipen=16)

rm(list=ls())


### parsing argument
args=commandArgs(trailingOnly=T)
if ( length(args)>0){
   trait <- args[1]
}else{
   trait <- "imputed_UKB_20002_1111_self_reported_asthma"
}   

### create directory 
outdir <- "./gwas_torusfile/"
if ( !file.exists(outdir) ) dir.create(outdir, showWarnings=F, recursive=T)

## traits <- sort(read.table("traits.txt")$V1)
 
######################################################
### Generate torus files with LD block information ###
######################################################


### read summary data
fn <- paste("../gwas_data/", trait, ".txt.gz", sep="")  
summ <- fread(fn, header=T, data.table=F)

if ( grepl("GBMI", trait)){
###
    zscore <- summ$inv_var_meta_beta/summ$inv_var_meta_sebeta
    summ <- cbind(summ[,c(1,2)], zscore)
    names(summ) <- c("chr", "pos", "zscore")
}else{
###    
   summ <- summ[,c(3,4,10)]
   names(summ) <- c("chr", "pos", "zscore")
}

###
summ <- summ%>%mutate(chr=as.integer(gsub("^chr", "", chr)), pos=as.integer(pos),
                      chr_pos=paste(chr, pos, sep="_"))%>%filter(chr%in%1:22)


### LD block files
bed <- read.table("eur_ld.hg38.bed", header=T) ##, fill=T, data.table=F)
bed <- bed%>%mutate(chr=as.integer(gsub("^chr", "", chr)))%>%filter(chr%in%1:22)%>%drop_na(start, stop)
 


### add locus using LD block file
nblock <- nrow(bed)
res_torus <- lapply(1:nblock, function(i){
   ##
   cat("Locus", i, "\n")
   ## 
   chr_i <- bed$chr[i]
   s0 <- as.numeric(bed$start[i])
   s1 <- as.numeric(bed$stop[i])

   ## SNPs in each block   
   summ2 <- summ%>%dplyr::filter(chr==chr_i, pos>=s0, pos<s1)
   if ( nrow(summ2)>0) { 
      summ2 <- summ2%>%mutate(Locus=paste("Loc", i, sep=""))%>%dplyr::select(id_b38=chr_pos, Locus, zscore)
   }else{
      summ2 <- NULL
   }
   summ2 
})
 
##
res_torus <- do.call(rbind, res_torus)

### output 
opfn <- paste0(outdir, trait, "_torus.txt.gz")
fwrite(res_torus, opfn, quote=F, sep=" ", row.names=F, col.names=F)


###END
