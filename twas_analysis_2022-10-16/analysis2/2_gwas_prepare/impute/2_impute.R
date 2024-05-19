###
library(data.table)
library(tidyverse)

options(scipen=16)

rm(list=ls())

### parsing argument
args=commandArgs(trailingOnly=T)
if ( length(args)>0 ){
   ##
   trait <- args[1]
   split_one <- args[2] 
}else{
   trait <- "Asthma_Bothsex_afr_inv_var_meta_GBMI_052021_nbbkgt1"
   split_one <- "splitGene000"
}       
    

###
### directory for output
outdir <- paste("./2_impute.outs/", trait, "/", sep="")
if ( !file.exists(outdir) ) dir.create(outdir, showWarnings=F, recursive=T)


## traits <- sort(read.table("traits.txt")$V1)

###
### gwas data, summary data
fn <- paste("../gwas_data/", trait, ".txt.gz", sep="")
res <- fread(fn, header=T, data.table=F)

##
if ( grepl("GBMI", trait)){
   ### 
   res2 <- res[,c(1:5, 7:9)]
   names(res2) <- c("chr", "pos", "ref", "alt", "rs", "meta_beta", "meta_se", "pval")    
   res2 <- res2%>%filter(chr%in%c(1:22))%>%
       mutate(id_variant=paste0("chr", paste(chr,  pos, ref, alt, "b38", sep="_")),
              chr_pos_grch38=paste(chr, pos, sep="_"), zscore=meta_beta/meta_se)
   res2 <- res2%>%dplyr::select(chr, pos, rs, id_variant, zscore, pval, chr_pos_grch38)  
}else{
   ###
   res2 <- res[,c(3:4, 1:2,  10:11)]
   names(res2) <- c("chr", "pos", "rs", "id_variant", "zscore", "pval")
   res2 <- res2%>%mutate(chr=as.integer(gsub("^chr", "", chr)), pos=as.integer(pos),
                         chr_pos_grch38=paste(chr, pos, sep="_"))%>%
       filter(chr%in%c(1:22))
}    
    
    
    

####
### missing snps
fn <- paste("./1_missing_snp.outs/", trait, "/", split_one, sep="")
snp <- read.table(fn)$V1
bed <- as.data.frame(str_split(snp, "_", simplify=T))
names(bed) <- c("chr", "pos")
bed$chr_pos_grch38 <- snp



###
### imputation 
nsnp <- nrow(bed)
time0 <- Sys.time()
tmp <- lapply(1:nsnp, function(i){
   ##
   if ( i%%10==0) cat(i, "\n")
   chr_i <- as.integer(bed$chr[i])
   pos_i <- as.integer(bed$pos[i])

   ### +/- 10 kb
   s0 <- pos_i-1e+04
   s1 <- pos_i+1e+04
    
   tmp2 <- res2%>%filter(chr==chr_i)
   pos <- tmp2$pos
   sub0 <- pos>=s0&pos<s1
   nsnp <- sum(sub0)

   if ( nsnp>0){
      ##
      tmp2 <- tmp2[sub0,]
      dtss <- abs(as.numeric(tmp2$pos)-pos_i)
      imin <- which.min(dtss)
      tmp3 <- tmp2[imin,]
      tmp3$id_b38_0 <- bed$chr_pos_grch38[i]
  }else{
      tmp3 <- NULL
  }    
  tmp3
})

###
###
time1 <- Sys.time()
diff0 <- difftime(time1, time0, units="secs")
cat(trait, diff0, "\n")

###
tmp <- do.call(rbind, tmp)

gfn <- paste(outdir, split_one, "_missing.txt.gz", sep="")
fwrite(tmp, file=gfn, sep="\t", col.names=F, quote=F, na=NA)

###
### output
## res3 <- res2%>%mutate(id_b38_0=chr_pos_grch38)
## res_comb <- rbind(res3, tmp)

## gfn <- paste(outdir, trait, "_impute.txt.gz", sep="")
## fwrite(res_comb, file=gfn, sep="\t", quote=F, na=NA)

###
### END









