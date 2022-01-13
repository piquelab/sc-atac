
###
library(tidyverse)
library(data.table)

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
### annotation files
anno <- read.table("/wsu/home/groups/piquelab/data/gencode/Gencode_human/release_31/GRCh37_mapping/gencode.v31lift37.annotation.gff3.gz", header=F, stringsAsFactors=F)

anno2 <- anno%>%
   mutate(gene_id=gsub("ID=|;.*", "", V9), ID=gsub("\\..*", "", gene_id))%>%
   dplyr::filter(V3=="gene")%>%
   dplyr::rename("Chr"="V1", "min"="V4", "max"="V5", "strand"="V7")%>%
   mutate(Chr=gsub("chr", "", Chr),
          start=ifelse(strand=="+", min, max),
          end=ifelse(strand=="-", max, min))

###
### vcf files
vcf <- fread("SCAIP1-6_infor.txt.gz", header=T, data.table=F, stringsAsFactors=F)

###
### response motif
snpanno <- fread("../SNPannotation/annot_jaspar2020_torus/Response_motif.annot.gz", header=T, data.table=F)


###
### prepare eQTL, gene map and snp map files for torus
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

### 2
### gene map file
   geneList <- unique(summ$gene)
   gene.map <- data.frame(gene=anno2$ID, chr=anno2$Chr, start=anno2$start, start2=anno2$start)%>%
      dplyr::filter(gene%in%geneList) 
   opfn <- gzfile(paste(outdir, condition, ".gene.map.gz", sep=""))
   write.table(gene.map, opfn, quote=F, row.names=F, col.names=F)
   ## close(opfn)

### 3   
### snp map file
   snpList <- unique(summ$SNP)
   snp.map <- data.frame(SNP=vcf$V3, chr=vcf$V1, pos=vcf$V2)%>%dplyr::filter(SNP%in%snpList)
   opfn <- gzfile(paste(outdir, condition, ".snp.map.gz", sep=""))
   write.table(snp.map, opfn, quote=F, row.names=F, col.names=F)
   ## close(opfn)

### 4
### annotation files
   snpanno2 <- snpanno%>%filter(SNP%in%res$V2)
   opfn <- gzfile(paste(outdir, condition, ".annot.gz", sep=""))
   write.table(snpanno2, opfn, quote=F, row.names=F, col.names=T)
   
}
###



###
### annotation one by one
### response motif

outdir <- "./torus_input/annot_positive/"
if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)

snpanno <- fread("../SNPannotation/annot_jaspar2020_torus/Positive_motif.annot.gz", header=T, data.table=F)
conditions <- colnames(snpanno)

###
for ( i in 2:ncol(snpanno)){
###
   anno2 <- snpanno[,c(1,i)]
   condition <- gsub("_d", "", conditions[i])  
   opfn <- gzfile(paste(outdir, "zzz_motif.", condition, ".annot.gz", sep=""))
   write.table(anno2, opfn, quote=F, row.names=F, col.names=T)
   cat(condition, "\n") 
}

## ###
## ### 
## summ <- summ%>%mutate(Chr=gsub(":.*","",SNP))
## summ2 <- summ%>%dplyr::filter(abs(t.stat)<2.5)
## opfn <- gzfile(paste(outdir, condition, ".filter.txt.gz", sep=""))
## write.table(summ2[,-6], opfn, quote=F, row.names=F) 


## datasets <- read.table("datasets_motif.txt")
## write.table(sort(datasets$V1), "datasets_motif.txt", quote=F, row.names=F, col.names=F)
