##
library(tidyverse)
library(data.table)
##
options(scipen=200)

### read summary data
fn <- "imputed_UKB_20002_1111_self_reported_asthma.txt.gz"
summ <- fread(fn, header=T, data.table=F)

summ <- summ%>%mutate(chr_pos=paste(gsub("chr", "", chromosome), position, sep="_"))%>%
   dplyr::select(variant_id, panel_variant_id, chr=chromosome,  pos=position, chr_pos, zscore)


### LD block files
bed <- fread("eur_ld.hg38.bed", header=T, fill=T, data.table=F)


## summ <- data.frame(SNP=res$V2, gene=res$V1, beta=res$V5, "t-stat"=zscore, "p-value"=res$V4)

###
nblock <- nrow(bed)
summ2 <- lapply(1:nblock, function(i){
   ##
   cat("Locus", i, "\n")
   ## 
   chr_i <- bed$chr[i]
   s0 <- as.numeric(bed$start[i])
   s1 <- as.numeric(bed$stop[i])
   ## 
   summ2 <- summ%>%dplyr::filter(chr==chr_i, as.numeric(pos)>=s0,as.numeric(pos)<s1)

   summ2 <- summ2%>%mutate(Locus=paste("Loc", i, sep="")) ##%>%dplyr::select(variant_id, Locus, zscore)
   summ2     
})
 
##
summ2 <- do.call(rbind, summ2)

### output 
opfn <- gzfile("Asthma_torus_zval.txt.gz", "w")
write.table(summ2, opfn, quote=F, row.names=F)

## x <- fread("Asthma_torus_zval.txt.gz", header=T, fill=T)


###
###
## id_panel <- as.character(summ[,2])
## names(id_panel) <- as.character(summ[,1])

## x <- fread("Asthma_torus_zval.txt.gz", header=T, data.table=F)

## x$SNP <- id_panel[as.character(x$SNP)]

## opfn <- gzfile("Asthma_torus_zval.txt.gz", "w")
## write.table(x, opfn, quote=F, sep="\t", row.names=F)
