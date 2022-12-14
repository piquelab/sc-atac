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


####################
### torus format ###
####################

### eQTL summary results
fn <- "./eQTL_results/Whole_Blood.allpairs.txt.gz"
res <- fread(fn, header=T, data.table=F, stringsAsFactors=F)
zscore <- res$slope/res$slope_se
summ <- data.frame(SNP=res$variant_id, gene=res$gene_id, beta=res$slope, "t-stat"=zscore, "p-value"=res$pval_nominal)
##
##summ2 <- summ%>%mutate(gene=gsub("\\..*", "", gene))

opfn <- gzfile(paste(outdir, "Whole_Blood.eQTL.txt.gz", sep=""))
write.table(summ, opfn, quote=F, row.names=F)





##########################
### gene map for torus ###
##########################

###
### annotation files
## anno <- read.table("/wsu/home/groups/piquelab/data/gencode/Gencode_human/release_31/GRCh37_mapping/gencode.v31lift37.annotation.gff3.gz", header=F, stringsAsFactors=F) ##grch37

###
### gene annotation files from grch38

fn <- "/wsu/home/groups/piquelab/data/gencode/Gencode_human/release_31/gencode.v31.annotation.gff3.gz"
## fn <- "/wsu/home/groups/piquelab/data/gencode/Gencode_human/release_42/gencode.v42.chr_patch_hapl_scaff.annotation.gff3.gz"
anno <- read.table(fn, header=F, stringsAsFactors=F)

anno2 <- anno%>%
   mutate(gene_id=gsub("ID=|;.*", "", V9), ID=gsub("\\..*", "", gene_id))%>%
   dplyr::filter(V3=="gene")%>%
   dplyr::rename("Chr"="V1", "min"="V4", "max"="V5", "strand"="V7")%>%
   mutate(Chr=gsub("chr", "", Chr),
          start=ifelse(strand=="+", min, max),
          end=ifelse(strand=="-", max, min))

###
### gene used for GTEx 
geneList <- read.table("./eQTL_results/geneList.txt")$V1 ### 20,315 gene
geneDf <- data.frame(gene_id2=unique(geneList))%>%mutate(ID=gsub("\\..*", "", gene_id2))

anno3 <- anno2%>%inner_join(geneDf, by="ID")%>%dplyr::filter(!grepl("PAR_Y", gene_id)) ## 20207


###
### gene map file
gene.map <- data.frame(gene=anno3$gene_id2, chr=anno3$Chr, start=anno3$start, start2=anno3$start)
###
opfn <- gzfile(paste(outdir, "zzz_gene.map.gz", sep=""))
write.table(gene.map, opfn, sep="\t", quote=F, row.names=F, col.names=F)




#########################################
### SNP map file and annotation files ###
#########################################

###
### GTEx SNP map file, grch38  
snpfile <- fread("./eQTL_results/snpList.txt", header=F, data.table=F) ### 9,719,090 million SNPs
snpfile2 <- snpfile%>%separate(V1, into=c("chr", "pos", "ref", "alt", "hg"), sep="_")%>%
   mutate(chr2=gsub("chr", "", chr)) 
snp.map <- data.frame(SNP=snpfile$V1, chr=snpfile2$chr2, pos=snpfile2$pos)## %>%
## snp.map <- snp.map%>%mutate(chr_pos=paste(chr, pos, sep="_"))

###
### our used geno file lift from grch37 to grch38
### V1-V4 grch37 version, V5 rs id including original ones and id recovered from annovar
### V6 grch38 version after liftover
bedlift38 <- fread("SCAIP_final_bed.gz", header=F, data.table=F)
DF_lift <- data.frame(id_b37=bedlift38$V4, chr_pos_grch38=bedlift38$V6)


snp.map2 <- snp.map%>%
    mutate(chr_pos_grch38=paste(chr, pos, sep="_"))%>%inner_join(DF_lift, by="chr_pos_grch38")
### 5,072,688 SNPs

id_37 <- snp.map2$id_b37
id_38 <- as.character(snp.map2$SNP)
names(id_38) <- as.character(id_37)

###
### combine2 SNP annotation files, multiple columns 
fn <- "/nfs/rprdata/julong/sc-atac/genetic.analysis_torus_2021-10-14/SNPannotation/4_SNPAnnot.outs/pct_0.02_combineNew/combine2_torus.annot.gz"

anno <- fread(fn, header=T, data.table=F)
anno2 <- anno%>%dplyr::filter(SNP%in%id_37)%>%mutate(SNP=id_38[as.character(SNP)]) ### 5,072,371

### output annotation
opfn <- gzfile(paste(outdir, "zzz_combine2_torus.annot.gz", sep=""))
write.table(anno2, opfn, quote=F, row.names=F, col.names=T)

## fn <- paste(outdir, "zzz_combine2_torus.annot.gz", sep="")
## x <- fread(fn, header=F)
## names(x) <- names(anno)
## gfn <- gzfile(paste(outdir, "zzz_combine2_torus.annot.gz", sep=""))
## write.table(x, gfn, quote=F, row.names=F, col.names=T)




###
### union annotation files
fn <- "/nfs/rprdata/julong/sc-atac/genetic.analysis_torus_2021-10-14/SNPannotation/4_SNPAnnot.outs/pct_0.02_combineNew/Union_torus.annot.gz"
anno <- fread(fn, header=T, data.table=F)
anno2 <- anno%>%dplyr::filter(SNP%in%id_37)%>%mutate(SNP=id_38[as.character(SNP)]) ### 5,072,371

### output annotation
opfn <- gzfile(paste(outdir, "zzz_union_torus.annot.gz", sep=""))
write.table(anno2, opfn, quote=F, row.names=F, col.names=T)

## fn <- paste(outdir, "zzz_union_torus.annot.gz", sep="")
## x <- fread(fn, header=F)
## names(x) <- c("SNP", "peaking_d")
## gfn <- gzfile(paste(outdir, "zzz_union_torus.annot.gz", sep=""))
## write.table(x, gfn, quote=F, row.names=F, col.names=T)


### output SNP map files
snp.map_final <- snp.map2%>%dplyr::filter(SNP%in%anno2$SNP)%>%dplyr::select(SNP, chr, pos) ## 5,072,433
opfn2 <- gzfile(paste(outdir, "zzz_snp.map.gz", sep=""))
write.table(snp.map_final, opfn2,  quote=F, row.names=F, col.names=F)


