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

 
bed <- fread("SCAIP_final_bed.gz", header=T, data.table=F)

id_snp <- bed$genetic_variant
names(id_snp) <- bed$chr_pos_grch38


###
traits <- read.table("./traits.txt")$V1[3:4]
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
fwrite(x3, file=opfn, row.names=F, col.names=F, sep=" ")    

###    
cat(ii, nrow(x), nrow(x3), "\n")
    
}    


###########################################
### if not tested in eqtl, discard it
############################################

outdir <- "./tmp/gwas_PIP2/"
if ( !file.exists(outdir) ) dir.create(outdir, showWarnings=F, recursive=T)


### replace SNP id and keep the same to those in eqtl mapping.

 
bed <- fread("SCAIP_final_bed.gz", header=T, data.table=F)

id_snp <- bed$genetic_variant
names(id_snp) <- bed$chr_pos_grch38


###
traits <- read.table("./traits.txt")$V1[3:4]
for (ii in traits){

###
cat(ii, "\n")
    
fn <- paste("./gwas_PIP/", ii, ".pip.gz", sep="")
x <- fread(fn, header=F, data.table=F)    
x <- x%>%distinct(V1, .keep_all=T)

x2 <- x%>%mutate(in_aloft=V1%in%bed$chr_pos_grch38,
                 id_aloft=ifelse(in_aloft, id_snp[V1], V1))%>%
    filter(in_aloft)
    
### output    
x3 <- x2%>%dplyr::select(id_aloft, V2, V3, V4)%>%drop_na(id_aloft, V2, V3, V4)    
opfn <- paste0(outdir, ii, ".pip.gz")
fwrite(x3, file=opfn, row.names=F, col.names=F, sep=" ")    

###    
cat(ii, nrow(x), nrow(x3), "\n")
    
}


 



    
## keep the same to old version    
## imputation gwas
## fn2 <- paste("./gwas_impute/", ii, "_impute.txt.gz", sep="")
## summ <- fread(fn2, header=T, data.table=F)    
## summ2 <- summ%>%dplyr::select(id_b38_0, chr_pos_grch38)

## ###    
## comb <- summ2%>%left_join(x, by=c("chr_pos_grch38"="V1"))
     
## comb2 <- comb%>%
##     mutate(in_aloft=id_b38_0%in%bed$chr_pos_grch38, id_aloft=ifelse(in_aloft, id_snp[id_b38_0], id_b38_0))%>%
##     filter(in_aloft)
    
## comb3 <- comb2%>%dplyr::select(id_aloft, V2, V3, V4)%>%drop_na(id_aloft, V2, V3, V4)
## opfn <- paste0(outdir, ii, "_old.pip.gz")    
## fwrite(comb3, file=opfn, row.names=F, col.names=F, sep=" ")



  
