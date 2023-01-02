###
library(tidyverse)
library(data.table)

rm(list=ls())

options(scipen=200)


### passing argument
## args=commandArgs(trailingOnly=T)
## if (length(args)>0){
##    condition <- args[1]
##  }else{
##    condition <- "Bcell_CTRL"
## }

outdir <- "./torus_input/"
if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)



################################
### asthma gwas summary data ###
################################

### add chr_pos information once and don't need run it again. 
fn <- "./gwas_input/Asthma_torus_zval.txt.gz"
summ <- fread(fn, header=T, data.table=F, fill=T)

###
summ2 <- summ%>%dplyr::select(SNP=panel_variant_id, Locus, zscore)
gfn <- gzfile("./torus_input/Asthma_all_torus_zval.txt.gz")
write.table(summ2, gfn, row.names=F, col.names=T, quote=F)


### chr, pos, chr_pos, variant id from GRCh37, rs id and chr_pos from GRCh38
## snp <- fread("SCAIP_final_bed.gz", header=F, data.table=F)

## ### based id
## ## share_SNP <- intersect(summ$SNP, snp$V5)
## ### length(share_SNP), 4,450,073

## ### based pos
## share_SNP <- intersect(summ$chr_pos, snp$V6)
## ### length(share_SNP), 4,632,548

## summ_final <- summ%>%filter(chr_pos%in%share_SNP)%>%
##    dplyr::select(-SNP, SNP=chr_pos, Locus, zscore)

## summ_final <- summ_final%>%dplyr::select(SNP, Locus, zscore)

## ###
## gfn <- gzfile("./torus_input/Asthma_torus_zval.txt.gz")   ### 4,632,554 SNPs
## write.table(summ_final, gfn, row.names=F, col.names=T, quote=F, sep="\t")

## summ <- fread("./torus_input/Asthma_torus_zval.txt.gz")

## based on rs, 4450079
##

### gwas
## summ2 <- summ%>%filter(SNP%in%share_SNPs)
## opfn2 <- gzfile("./torus_input/Asthma_torus_zval.txt.gz", "w")
## write.table(summ2, opfn2, row.names=F, quote=F, sep="\t", col.names=T)
## close(opfn2)



###
### annotation files
## MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
## for (i in 1:3){
## ###    
## oneMCl <- MCls[i]

## cat(oneMCl, "\n")
    
## fn <- paste("../SNPannotation/4_SNPAnnot.outs/pct_0.02/3_", oneMCl, "_union_torus.annot.gz", sep="")
## annot <- fread(fn, header=T, data.table=F)
## annot <- annot%>%mutate(rs=gsub(".*;", "", SNP))

## share_SNPs <- intersect(summ$SNP, annot$rs)

## ###
## ### output

## ## annotation
## annot2 <- annot%>%filter(rs%in%share_SNPs)%>%dplyr::select(SNP=rs, peaking_d)
## opfn <- gzfile(paste("./torus_input/", oneMCl, "_union_torus_annot.txt.gz", sep=""), "w")
## write.table(annot2, opfn, row.names=F, quote=F, sep="\t", col.names=T)
## close(opfn)

## }

####################################################################
### annotation files, only annotate genetic variants within peak ###
####################################################################

outdir2 <- "./torus_input/Peak_annot/"
if ( !file.exists(outdir2)) dir.create(outdir2, showWarnings=F, recursive=T)


summ <- fread("./torus_input/Asthma_torus_zval.txt.gz")
bed <- fread("SCAIP_final_bed.gz", header=F, data.table=F)
bed <- bed%>%filter(V6%in%summ$SNP)

chr_pos <- bed$V6
names(chr_pos) <- as.character(bed$V4)


MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
for (i in 1:4){
###    
oneMCl <- MCls[i]

cat(oneMCl, "\n")
    
fn <- paste("../SNPannotation/4_SNPAnnot.outs/pct_0.02/", oneMCl, "_Active_peak_torus.annot.gz", sep="")
annot <- fread(fn, header=T, data.table=F)
annot2 <- annot%>%filter(SNP%in%bed$V4)%>%mutate(chr_pos_grch38=chr_pos[as.character(SNP)])
annot2 <- annot2%>%dplyr::select(chr_pos_grch38, peaking_d)%>%dplyr::rename(SNP=chr_pos_grch38)

## output
gfn <- gzfile(paste(outdir2, oneMCl, "_peaks_torus.annot.gz", sep=""))
write.table(annot2, gfn, row.names=F, quote=F, sep="\t", col.names=T)
##
}



##########################
### combine annotation ###
##########################


outdir2 <- "./torus_input/Combine/"
if ( !file.exists(outdir2)) dir.create(outdir2, showWarnings=F, recursive=T)


## summ <- fread("./torus_input/Asthma_torus_zval.txt.gz")
bed <- fread("SCAIP_final_bed.gz", header=F, data.table=F)
## bed <- bed%>%filter(V6%in%summ$SNP)
DF_lift <- data.frame(id_b37=as.character(bed$V4), chr_pos_grch38=as.character(bed$V6))

summ <- fread("./gwas_input/Asthma_torus_zval.txt.gz", header=T, data.table=T)

summ2 <- summ%>%inner_join(DF_lift, by=c("chr_pos"="chr_pos_grch38"))

## x <- summ2[1:20, c("panel_variant_id", "chr_pos", "id_b37")]
## openxlsx::write.xlsx(x, "toy_asthma.xlsx", overwrite=T)


##bed <- fread("SCAIP_final_bed.gz", header=F, data.table=F)
chr_pos_grch38 <- as.character(bed$V6)
names(chr_pos_grch38) <- as.character(bed$V4)


###
### multiple column annotation 
fn <- "/nfs/rprdata/julong/sc-atac/genetic.analysis_torus_2021-10-14/SNPannotation/4_SNPAnnot.outs/pct_0.02_combineNew/combine2_torus.annot.gz"
anno <- fread(fn, header=T, data.table=F)

anno2 <- anno%>%filter(SNP%in%as.character(bed$V4))%>%
    mutate(SNP=chr_pos_grch38[as.character(SNP)])

gfn <- gzfile(paste(outdir2, "combine2_torus.annot.gz", sep=""))
write.table(anno2, file=gfn, quote=F, row.names=F)



###
###
## fn <- "./torus_input/Combine/combine2_torus.annot.gz"
## anno <- fread(fn, header=T, data.table=F)


###
### single column annotation, peaks, hit by cell-type motifs and hit by response motifs
fn <- "/nfs/rprdata/julong/sc-atac/genetic.analysis_torus_2021-10-14/SNPannotation/4_SNPAnnot.outs/pct_0.02_combineNew/Union_torus.annot.gz"
anno <- fread(fn, header=T, data.table=F)

anno2 <- anno%>%filter(SNP%in%as.character(bed$V4))%>%
    mutate(SNP=chr_pos_grch38[as.character(SNP)])

gfn <- gzfile(paste(outdir2, "Union_torus.annot.gz", sep=""))
write.table(anno2, file=gfn, quote=F, row.names=F)



    


## End, Oct-2-2022











       
 
