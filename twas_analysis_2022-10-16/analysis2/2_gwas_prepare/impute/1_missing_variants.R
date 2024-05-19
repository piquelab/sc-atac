##
library(tidyverse)
library(data.table)
library(openxlsx)
##
options(scipen=16)

rm(list=ls())

outdir <- "./1_missing_snp.outs/"
if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)


###
### aloft variants directory, include grch38 position 
fn <- "../SCAIP_final_bed.gz"
snp_df <- fread(fn, header=T, data.table=F)



###
### ALOFT missing
fn <- "/nfs/rprdata/julong/sc-atac/twas_analysis_2022-10-16/ALOFT_results/aloft_grch37.txt.gz"
topsnp <- read.table(fn, header=F)$V1
df_aloft <- snp_df%>%filter(genetic_variant%in%topsnp)


###
### GTEx
fn2 <- "/nfs/rprdata/julong/sc-atac/twas_analysis_2022-10-16/Whole_Blood_GTEx_v8_results/gtex_grch38.txt.gz"
topsnp2 <- read.table(fn2, header=F)$V1

x <- str_split(gsub("^chr", "", topsnp2), "_", simplify=T)    
df_gtex <- data.frame(chr=x[,1], chr_pos=paste(x[,1], x[,2], sep="_"))%>%
    filter(chr%in%as.character(1:22))

  
traits <- sort(read.table("traits.txt")$V1)[3:4]
summary_df <- NULL
for (ii in traits){

###
### summary data
   fn <- paste("../gwas_data/", ii, ".txt.gz", sep="")
   summ <- fread(fn, header=T, data.table=F)

   if ( grepl("GBMI", ii)){
   ### 
      x <- summ[,1:2]
      names(x) <- c("chr", "pos")
      x2 <- x%>%mutate(chr_pos=paste(chr, pos, sep="_"))%>%filter(chr%in%as.integer(1:22))
      snp_gwas <- unique(x2$chr_pos)
       
   }else{
   ###
      x <- str_split(gsub("^chr", "", summ[,2]), "_", simplify=T)
      x2 <- data.frame(chr=x[,1], chr_pos=paste(x[,1], x[,2], sep="_"))%>%filter(chr%in%as.character(1:22))
      snp_gwas <- unique(x2$chr_pos)    
   }    

   snp_miss_aloft <- df_aloft%>%filter(!chr_pos_grch38%in%snp_gwas)%>%pull(chr_pos_grch38)%>%unique()
   snp_miss_gtex <- df_gtex%>%filter(!chr_pos%in%snp_gwas)%>%pull(chr_pos)%>%unique()
   snp_miss <- union(snp_miss_aloft, snp_miss_gtex)   
        
   ###
   miss_aloft <- length(snp_miss_aloft)
   n_aloft <- nrow(df_aloft)
   miss_gtex <- length(snp_miss_gtex)
   n_gtex <- nrow(df_gtex)
   miss <- length(snp_miss)
    
   tmp0 <- data.frame(traits=gsub("_inv_.*", "", ii), nsnp_gwas=nrow(summ),
                      total_miss=miss,  
                      nsnp_aloft=paste0(miss_aloft, "|", n_aloft),
                      nsnp_gtex=paste0(miss_gtex, "|", n_gtex))
  summary_df <- rbind(summary_df, tmp0)  
  cat(ii, length(snp_miss), "\n")    
    
  ### output
  tmp <- str_split(snp_miss, "_", simplify=T)
  if ( ncol(tmp)>2){
  ### 
     isel <- sapply(3:ncol(tmp), function(k) tmp[,k]=="")
     isel2 <- apply(isel, 1, all)
     ###
     snp_miss <- snp_miss[isel2]    
  }
   
  opfn <- paste(outdir, ii,  "_missing.txt", sep="")
  write.table(snp_miss, file=opfn, quote=F, row.names=F, col.names=F)       
}    


### output
opfn <- paste(outdir, "summary_snp.txt", sep="")
write.table(summary_df, file=opfn, row.names=F, col.names=T, sep="\t")


###
### End


##############################
### summary missing SNPs 
#############################



traits <- sort(read.table("traits.txt")$V1)

###
### aloft variants directory, include grch38 position 
fn <- "../SCAIP_final_bed.gz"
bed <- fread(fn, header=T, data.table=F)

## snp_df <- snp_df%>%mutate(chr_grch38=gsub("_.*", "", chr_pos_grch38))

ii <- traits[4]
fn <- paste("../gwas_data/", ii, ".txt.gz", sep="")
summ <- fread(fn, header=T, data.table=F)

if ( grepl("GBMI", ii)){
### 
   x <- summ[,1:2]
   names(x) <- c("chr", "pos")
   x2 <- x%>%mutate(chr_pos=paste(chr, pos, sep="_"))%>%filter(chr%in%as.integer(1:22))
   snp_gwas <- unique(x2$chr_pos)
       
}else{
###
   x <- str_split(gsub("^chr", "", summ[,2]), "_", simplify=T)
   x2 <- data.frame(chr=x[,1], chr_pos=paste(x[,1], x[,2], sep="_"))%>%filter(chr%in%as.character(1:22))
   snp_gwas <- unique(x2$chr_pos)    
}    

old <- unique(bed$chr_pos_grch38)
olap <- intersect(old, snp_gwas)



## ### ALOFT
## topsnp <- read.table("aloft_grch37.txt.gz", header=F)$V1

## ### GTEx
## topsnp2 <- read.table("gtex_grch38.txt.gz", header=F)$V1

 

## ###
## summ <- map_dfr(traits, function(ii){
## ###
## ### summary data
## fn <- paste("../", ii, ".txt.gz", sep="")
## res <- fread(fn, header=T, data.table=F)

## ##
## res2 <- res[,c(1:5,7:9)]
## names(res2) <- c("chr", "pos", "ref", "alt", "rs", "meta_beta", "meta_se", "meta_p")
## res2 <- res2%>%filter(chr%in%c(1:22))%>%mutate(chr_pos_grch38=paste(chr, pos, sep="_"))
## res3 <- res2%>%drop_na(rs)

## ###
## ### ALOFT missing 
## x <- str_split(gsub(";.*", "", topsnp), ":", simplify=T)
## snp2 <- paste(x[,1], x[,2], sep="_")
## snp_miss <- snp_df%>%
##     filter(chr_pos%in%snp2, !chr_pos_grch38%in%res3$chr_pos_grch38, chr_grch38%in%as.character(1:22))%>%
##     pull(chr_pos_grch38)%>%unique()

## tmp <- str_split(snp_miss, "_", simplify=T)
## if ( ncol(tmp)>2){
##    ### 
##    isel <- sapply(3:ncol(tmp), function(k) tmp[,k]=="")
##    isel2 <- apply(isel, 1, all)
##    ###
##    snp_miss <- snp_miss[isel2]    
## }
    


## ###
## ### GTEx
## ## topsnp2 <- read.table("gtex_grch38.txt.gz", header=F)$V1
## x <- str_split(gsub("^chr", "", topsnp2), "_", simplify=T)    
## x2 <- data.frame(chr=x[,1], chr_pos=paste(x[,1], x[,2], sep="_"))
## snp_miss2 <- x2%>%filter(chr%in%as.character(1:22), !chr_pos%in%res3$chr_pos_grch38)%>%pull(chr_pos)%>%unique() 

## tmp <- str_split(snp_miss2, "_", simplify=T)
## if ( ncol(tmp)>2){
##    ### 
##    isel <- sapply(3:ncol(tmp), function(k) tmp[,k]=="")
##    isel2 <- apply(isel, 1, all)
##    ###
##    snp_miss2 <- snp_miss2[isel2]    
## }

## n_aloft <- length(topsnp)
## miss_aloft <- length(snp_miss)
## n_gtex <- length(topsnp2)
## miss_gtex <- length(snp_miss2)    
## summ0 <- data.frame(traits=gsub("_inv_.*", "", ii), nsnp_gwas=nrow(res3),
##                     nsnp_aloft=paste0(miss_aloft, "|", n_aloft),
##                     nsnp_gtex=paste0(miss_gtex, "|", n_gtex))
## cat(ii, "\n")
## summ0       
## })    

 
## ### output
## opfn <- paste(outdir, "summary_snp.xlsx", sep="")
## write.xlsx(summ, file=opfn)
