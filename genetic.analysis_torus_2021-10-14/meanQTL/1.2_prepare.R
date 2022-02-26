#
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
