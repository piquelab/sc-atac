### script prepare input files for torus
### create by JW, 2/21/2022

library(tidyverse)
library(data.table)

rm(list=ls())

### passing argument
## args=commandArgs(trailingOnly=T)
## if (length(args)>0){
##    condition <- args[1]
##  }else{
##    condition <- "Bcell_CTRL"
## }

outdir <- "./torus_input/"
if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)



###
### prepare eQTL for torus
datasets <- read.table("datasets.txt", header=F)
for (condition in datasets$V1){
   cat(condition, "\n")
### 1   
### eQTL summary results   
   fn <- paste("./eQTL_results/", condition, ".nominals.all.chunks.txt.gz", sep="")
   res <- fread(fn, sep=" ", data.table=F, stringsAsFactors=F)
   zscore <- abs(qnorm(res$V4/2))*sign(res$V5)
   summ <- data.frame(SNP=res$V2, gene=res$V1, beta=res$V5, "t-stat"=zscore, "p-value"=res$V4)
   ## summ <- summ%>%dplyr::filter("t-stat"!=0)
   opfn <- gzfile(paste(outdir, condition, ".eQTL.txt.gz", sep=""))
   write.table(summ, opfn, quote=F, row.names=F)
   ## close(opfn)
}
###


## ### 2
## ### gene map file
##    geneList <- unique(summ$gene)
##    gene.map <- data.frame(gene=anno2$ID, chr=anno2$Chr, start=anno2$start, start2=anno2$start)%>%
##       dplyr::filter(gene%in%geneList) 
##    opfn <- gzfile(paste(outdir, condition, ".gene.map.gz", sep=""))
##    write.table(gene.map, opfn, quote=F, row.names=F, col.names=F)
##    ## close(opfn)

## ### 3   
## ### snp map file
##    snpList <- unique(summ$SNP)
##    snp.map <- data.frame(SNP=vcf$V3, chr=vcf$V1, pos=vcf$V2)%>%dplyr::filter(SNP%in%snpList)
##    opfn <- gzfile(paste(outdir, condition, ".snp.map.gz", sep=""))
##    write.table(snp.map, opfn, quote=F, row.names=F, col.names=F)
##    ## close(opfn)

## ### 4
## ### annotation files
##    snpanno2 <- snpanno%>%filter(SNP%in%res$V2)
##    opfn <- gzfile(paste(outdir, condition, ".annot.gz", sep=""))
##    write.table(snpanno2, opfn, quote=F, row.names=F, col.names=T)
   
## }
###


############################################
### gene map and SNP map files for torus ###
############################################

###
### annotation files
anno <- read.table("/wsu/home/groups/piquelab/data/gencode/Gencode_human/release_31/GRCh37_mapping/gencode.v31lift37.annotation.gff3.gz", header=F, stringsAsFactors=F)

### tested gene list
load("/nfs/rprdata/julong/SCAIP/analyses/SCAIP-B1-6_2020.03.23/10_RNA.Variance_output/tmp9/1.2_Sel.PhxNew.RData")
geneList <- gsub("\\..*", "", rownames(PhxNew2))
write.table(geneList, "zzz_geneList.txt", row.names=F, col.names=F, quote=F)

anno2 <- anno%>%
   mutate(gene_id=gsub("ID=|;.*", "", V9), ID=gsub("\\..*", "", gene_id))%>%
   dplyr::filter(V3=="gene")%>%
   dplyr::rename("Chr"="V1", "min"="V4", "max"="V5", "strand"="V7")%>%
   mutate(Chr=gsub("chr", "", Chr),
          start=ifelse(strand=="+", min, max),
          end=ifelse(strand=="-", max, min))

anno2 <- anno2%>%dplyr::filter(ID%in%geneList)


### vcf files
vcf <- fread("SCAIP1-6_infor.txt.gz", header=F, data.table=F, stringsAsFactors=F)


### 1
### SNP list
conditions <- read.table("datasets.txt")$V1
SNPlist <- NULL
for (ii in conditions){
    ###
    cat(ii, "\n") 
    fn <- paste("./eQTL_results/", ii, ".nominals.all.chunks.txt.gz", sep="")
    res <- fread(fn, sep=" ", data.table=F, stringsAsFactors=F)
    ##
    SNPs <- unique(res$V2)
    SNPlist <- unique(c(SNPlist, SNPs))
}
##
write.table(SNPlist,"./eQTL_results/zzz_SNPlist.txt", row.names=F, col.names=F, quote=F)


### 2
### gene map file
gene.map <- data.frame(gene=anno2$ID, chr=anno2$Chr, start=anno2$start, start2=anno2$start)%>%
   dplyr::filter(gene%in%geneList) 
opfn <- gzfile(paste(outdir, "zzz_gene.map.gz", sep=""))
write.table(gene.map, opfn, sep="\t", quote=F, row.names=F, col.names=F)


### 3   
### snp map file
snp.map <- data.frame(SNP=vcf$V3, chr=vcf$V1, pos=vcf$V2)## %>%
   ## dplyr::filter(SNP%in%SNPlist)%>%distinct(SNP,.keep_all=T)
###
opfn <- gzfile(paste(outdir, "zzz_snp.map.gz", sep=""))
write.table(snp.map, opfn, sep="\t", quote=F, row.names=F, col.names=F)




######################
### fastqtl format ###
######################


## outdir <- "./torus_input/fastq_format/"
## if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)

## datasets <- read.table("datasets.txt", header=F)
## for (condition in datasets$V1){
##    cat(condition, "\n")
## ### 1   
## ### eQTL summary results   
##    fn <- paste("./eQTL_results/", condition, ".nominals.all.chunks.txt.gz", sep="")
##    res <- fread(fn, sep=" ", data.table=F, stringsAsFactors=F)
##    zscore <- abs(qnorm(res$V4/2))*sign(res$V5)
##    sdbeta <- abs(res$V5/zscore)
##    res$V6 <- sdbeta
##    ## summ <- data.frame(gene=res$V1, SNP=res$V2, dtss=res$V3, beta=res$V5, "t-stat"=zscore, "p-value"=res$V4)
##    ## summ <- summ%>%dplyr::filter("t-stat"!=0)
##    opfn <- gzfile(paste(outdir, condition, ".eQTL.txt.gz", sep=""))
##    write.table(res, opfn, sep="\t", quote=F, row.names=F, col.names=T)
## }




