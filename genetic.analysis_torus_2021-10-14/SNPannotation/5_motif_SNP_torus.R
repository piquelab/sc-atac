##
library(tidyverse)
library(data.table)
## library(Seurat)
## library(Signac)
library(Matrix)

rm(list=ls())

## passing argument
args <- commandArgs(trailingOnly=T)
if ( length(args)>0){
    motifFile <- args[1]
   oneMCl <- args[2] 
}else{
    motifFile <- "splitMotif000"
    oneMCl <- "Tcell"
}    


####################
### main program ###
####################

fn <- "SCAIP1-6_infor.txt.gz"
allsnp <- fread(fn, header=F, data.table=F)%>%mutate(chr_pos=paste(V1, V2, sep="_"))

    
outdir <- paste("./4_SNPAnnot.outs/pct_0.02_", oneMCl, "/",  sep="")
if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)

    
## cell-type active peaks annotation    
fn <- paste("./4_SNPAnnot.outs/pct_0.02/",  oneMCl, "_Active_peak_torus.annot.gz", sep="")
anno <- fread(fn, header=T, data.table=F)


### read motif 
fn <- paste("./Motif_file/", motifFile, sep="")
motifList <- read.table(fn)$V1

### 
for (i in 1:length(motifList)){

    motif <- motifList[i]

   ### motif annotation file   
   fn <- paste("./annot_jaspar2020/allsnp_", motif, ".bed.gz", sep="")
   snp.motif <- try(fread(fn, header=F, data.table=F, stringsAsFactors=F), silent=T)


### 
if ( class(snp.motif)!="try-error"){

   ###
   snp.motif <- snp.motif%>%mutate(chr_pos=paste(V1, V2, sep="_"))
   snpSel <- allsnp%>%dplyr::filter(chr_pos%in%snp.motif$chr_pos)%>%dplyr::pull(V3)
   ###
   anno[anno[,2]==1&anno[,1]%in%snpSel,2] <- 2

   len <- sum(anno[,2]==2)
   if ( len>0){
   ###
      gfn <- gzfile(paste(outdir, motif, "_union_torus.annot.gz", sep=""))
      write.table(anno, gfn, sep="\t", row.names=F, col.names=T, quote=F)
   ###    
   }    
} ##    

}
### End
