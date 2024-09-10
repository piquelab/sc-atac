###
###
library(tidyverse)
library(data.table)

options(scipen=16)


######################################################################
### if not tested in eqtl, keep it with new id, default input data 
#######################################################################

outdir <- "./gwas_PIP2/"
if ( !file.exists(outdir) ) dir.create(outdir, showWarnings=F, recursive=T)


### replace SNP id and keep the same to those in eqtl mapping.

 
bed <- fread("gtex_v8_snpinfor.txt.gz", header=F, data.table=F)
bed <- bed[,1:3]
names(bed) <- c("chr", "pos", "genetic_variant")
bed <- bed%>%mutate(chr_pos_grch38=paste(gsub("^chr", "", chr), pos, sep="_"))

id_snp <- bed$genetic_variant
names(id_snp) <- bed$chr_pos_grch38


###
traits <- sort(read.table("./traits.txt")$V1)[3]
for (ii in traits){

###
cat(ii, "\n")
    
fn <- paste("./gwas_PIP/", ii, ".pip.gz", sep="")
x <- fread(fn, header=F, data.table=F)    
x <- x%>%distinct(V1, .keep_all=T)

x2 <- x%>%mutate(in_aloft=V1%in%bed$chr_pos_grch38,
                 id_aloft=ifelse(in_aloft, id_snp[V1], V1))
### output    
x3 <- x2%>%dplyr::select(id_aloft, V2, V3, V4)%>%drop_na(id_aloft, V2, V3, V4)    
opfn <- paste0(outdir, ii, ".pip.gz")
fwrite(x3, file=opfn, row.names=F, col.names=F, sep="\t")    

###    
cat(ii, nrow(x), nrow(x3), "\n")
    
}    



  
