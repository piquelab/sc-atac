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



### TORUS format  
### eQTL summary results   
fn <- "./eQTL_results/PC1-18.nominals.eQTL.txt.gz"
res <- fread(fn, sep=" ", data.table=F, stringsAsFactors=F)
zscore <- abs(qnorm(res$V4/2))*sign(res$V5)
summ <- data.frame(SNP=res$V2, gene=res$V1, beta=res$V5, "t-stat"=zscore, "p-value"=res$V4)
## summ <- summ%>%dplyr::filter("t-stat"!=0)
opfn <- gzfile(paste(outdir, "ALOFT_eQTL.txt.gz", sep=""))
write.table(summ, opfn, quote=F, row.names=F, sep="\t")


geneList <- unique(res$V1)

## opfn <- "zzz_geneList.txt"
## write.table(geneList, opfn, row.names=F, col.names=F, quote=F)



###
### annotation files
anno <- read.table("/wsu/home/groups/piquelab/data/gencode/Gencode_human/release_31/GRCh37_mapping/gencode.v31lift37.annotation.gff3.gz", header=F, stringsAsFactors=F)

### tested gene list
## load("/nfs/rprdata/julong/SCAIP/analyses/SCAIP-B1-6_2020.03.23/6_DEG.CelltypeNew_output/Filter2/YtX_sel.comb.RData")
## geneList <- gsub("\\..*", "", rownames(YtX_sel))
 
anno2 <- anno%>%
   mutate(gene_id=gsub("ID=|;.*", "", V9), ID=gsub("\\..*", "", gene_id))%>%
   dplyr::filter(V3=="gene")%>%
   dplyr::rename("Chr"="V1", "min"="V4", "max"="V5", "strand"="V7")%>%
   mutate(Chr=gsub("chr", "", Chr),
          start=ifelse(strand=="+", min, max),
          end=ifelse(strand=="-", max, min))

anno2 <- anno2%>%dplyr::filter(ID%in%geneList)



###
### vcf files
vcf <- fread("ALOFT_infor.txt.gz", header=F, data.table=F, stringsAsFactors=F)


###
### prepare eQTL, gene map and snp map files for torus
## datasets <- read.table("datasets.txt", header=F)$V1
## ## datasets <- c("Bcell", "Monocyte", "NKcell", "Tcell")
## indir <- "./eQTL_results/"
## for (condition in datasets){
##    cat(condition, "\n")

## ### 1   
## ### eQTL summary results   
##    fn <- paste(indir, condition, ".nominals.all.chunks.txt.gz", sep="")
##    res <- fread(fn, sep=" ", data.table=F, stringsAsFactors=F)
##    zscore <- abs(qnorm(res$V4/2))*sign(res$V5)
##    summ <- data.frame(SNP=res$V2, gene=res$V1, beta=res$V5, "t-stat"=zscore, "p-value"=res$V4)
##    ## summ <- summ%>%dplyr::filter("t-stat"!=0)
##    opfn <- gzfile(paste(outdir, condition, ".eQTL.txt.gz", sep=""))
##    write.table(summ, opfn, quote=F, row.names=F, sep="\t")
##    ## close(opfn)
## }


### 2
### gene map file
gene.map <- data.frame(gene=anno2$ID, chr=anno2$Chr, start=anno2$start, start2=anno2$start)
opfn <- gzfile(paste(outdir, "zzz_gene.map.gz", sep=""))
write.table(gene.map, opfn, quote=F, row.names=F, col.names=F, sep="\t")
   ## close(opfn)

### 3   
### snp map file
snp.map <- data.frame(SNP=vcf$V3, chr=vcf$V1, pos=vcf$V2)
opfn <- gzfile(paste(outdir, "zzz_snp.map.gz", sep=""))
write.table(snp.map, opfn, quote=F, row.names=F, col.names=F, sep="\t")
   ## close(opfn)


##########################
### summary annotation ###
##########################

conditions <- read.table("datasets.txt", header=F)$V1
## datasets <- c("Bcell", "Monocyte", "NKcell", "Tcell")
indir <- "./eQTL_results/"
summ2 <- map_dfr(1:20, function(i){
   ### 
   condition <- conditions[i]
   cat(i, condition, "\n") 
   ### eQTL summary results   
   fn <- paste(indir, condition, ".nominals.all.chunks.txt.gz", sep="")
   res <- fread(fn, sep=" ", data.table=F, stringsAsFactors=F)
   gene <- unique(res$V1)
   snps <- unique(res$V2)

   ## annotation
   oneMCl <- gsub("_.*", "", condition) 
   fn2 <- paste("../SNPannotation/4_SNPAnnot.outs/pct_0.02/3_", oneMCl, "_union_torus.annot.gz", sep="")
   annot <- fread(fn2, header=T, data.table=F)
   annot2 <- annot%>%filter(SNP%in%snps)
   ###
   df2 <- data.frame(conditions=condition, ngenes=length(gene),
        nsnp1=sum(annot2[,2]==1),
        nsnp2=sum(annot2[,2]==2),
        nsnp3=sum(annot2[,2]==3))
   df2
})

opfn <- "./torus_input/zzz_summary.snps.xlsx"
openxlsx::write.xlsx(summ2, opfn, overrite=T)

