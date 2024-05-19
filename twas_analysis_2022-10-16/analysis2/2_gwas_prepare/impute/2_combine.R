###
library(data.table)
library(tidyverse)

options(scipen=16)

rm(list=ls())

    

###
### directory for output
outdir <- "./2_impute.outs/"
if ( !file.exists(outdir) ) dir.create(outdir, showWarnings=F, recursive=T)


traits <- sort(read.table("traits.txt")$V1)[3:4]


###
for ( trait in traits){

cat(trait, "\n")    
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


res3 <- res2%>%mutate(id_b38_0=chr_pos_grch38)


fn2 <- paste(outdir, trait, "_all_missing.txt.gz", sep="")
tmp <- fread(fn2, header=F, data.table=F)
names(tmp) <- names(res3)
    
### output 
res_comb <- rbind(res3, tmp)
gfn <- paste(outdir, trait, "_impute.txt.gz", sep="")
fwrite(res_comb, file=gfn, sep="\t", quote=F, na=NA)

}


###
### END
