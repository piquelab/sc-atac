###
###
library(tidyverse)
library(data.table)
## library(Seurat)
## library(Signac)
library(Matrix)
 
### annotation SNP hit by peak, cell-type motif and response motif
### edit by Julong wei, 2-8-2022



### 1
### generate cell-type active peak bed file

fn <- "/nfs/rprdata/julong/sc-atac/analyses.2021-02-05/3_motif/1.3_motif.outs/1_cell-type_active.peaks.rds"
peak <- read_rds(fn)

MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
for(i in 1:length(MCls)){
   ###
   oneMCl <- MCls[i] 
   cat(oneMCl, "\n")
    
   outdir <- paste("./3_SNPAnnot.outs/", oneMCl,"/",  sep="")
   if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)
   ###
    
   peak2 <- peak[peak[,i+1]==1,] 
   peakSel <- peak2$peakAll
    
   bed <- str_split(peakSel, "-", simplify=T)

   opfn <- paste(outdir, oneMCl, "_Active_peak.bed", sep="")
   write.table(bed, opfn, sep="\t", row.names=F, quote=F, col.names=F)
}


### 2
### test if SNP is hit by cell-type active peak
## 3.2_intersect.sh


### 3
### torus format 
fn <- "SCAIP1-6_infor.txt.gz"
allsnp <- fread(fn, header=F, data.table=F)%>%mutate(chr_pos=paste(V1, V2, sep=":"))

MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
for( oneMCl in MCls){
###
   cat(oneMCl, "\n")
   outdir <- paste("3_SNPAnnot.outs/", oneMCl, "/", sep="")
   fn <- paste(outdir, oneMCl, "_Active_peak_snp.bed.gz", sep="")
   snpPeak <- fread(fn, header=F, data.table=F)
   snpPeak2 <- snpPeak%>%mutate(chr_pos=paste(V1, V2, sep=":"))%>%dplyr::filter(V4>0)%>%dplyr::pull(chr_pos) 
   peaking_d <- ifelse(allsnp$chr_pos%in%snpPeak2, 1, 0)
   anno <- data.frame("SNP"=allsnp$V3, "peaking_d"=peaking_d)
   ##
   opfn <- gzfile(paste(outdir, oneMCl, "_Active_peak_torus.annot.gz", sep=""))
   write.table(anno, opfn, sep="\t", row.names=F, col.names=T, quote=F) 
}


###
###
## testsnp <- read.table("../eQTL/eQTL_results/zzz_SNPlist.txt")$V1
## MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
## summ <- map_dfr(MCls, function(ii){
##     ##
##     fn <- paste(outdir, ii, "_peak.bed", sep="")
##     peak <- read.table(fn)
##     ##
##     fn <- paste(outdir, ii, "_torus_peak.annot.gz", sep="")
##     anno <- fread(fn, header=T, data.table=F)
##     anno <- anno%>%filter(SNP%in%testsnp)
##     ##
##     summ0 <- data.frame(MCls=ii, npeak=nrow(peak), nsnp=sum(anno[,2]))
##     summ0
## })    


### 4
### add motif and peak annotation ###

rm(list=ls())


fn <- "SCAIP1-6_infor.txt.gz"
allsnp <- fread(fn, header=F, data.table=F)%>%mutate(chr_pos=paste(V1, V2, sep=":"))

###
testsnp <- read.table("../eQTL/eQTL_results/zzz_SNPlist.txt")$V1


###
### annotation
oneMCl <- "Tcell"
fn <-"./3_SNPAnnot.outs/Tcell/Tcell_Active_peak_torus.annot.gz"
anno <- fread(fn, header=T, data.table=F)
anno <- anno%>%filter(SNP%in%testsnp)


###
### cell-type specific motif 
fn <- "/nfs/rprdata/julong/sc-atac/analyses.2021-02-05/3_motif/1.3_motif.outs/2_4_Tcell.selectMotif.txt" 
motifList <- read.table(fn, header=T)[,1]
snpmotif <- lapply(motifList, function(motif){
###
## cat(motif, "\n") 
   fn <- paste("./annot_jaspar2020/allsnp_", motif, ".bed.gz", sep="")   
   snp.motif <- try(fread(fn, header=F, data.table=F, stringsAsFactors=F), silent=T)
   if ( class(snp.motif)!="try-error"){ 
      snp.motif <- snp.motif%>%mutate(chr_pos=paste(V1, V2, sep=":"))
      xx <- snp.motif$chr_pos
   }else{
      xx <- NA
   }  
   xx  
})
snpmotif <- unlist(snpmotif[!is.na(snpmotif)])
snpmotif <- unique(snpmotif)
snpmotif2 <- allsnp%>%filter(chr_pos%in%snpmotif)%>%dplyr::pull(V3)
snpmotif2 <- unique(snpmotif2)

anno[anno[,2]==1&anno[,1]%in%snpmotif2, 2] <- 2 


###
### response motifs
outdir <- "./3_SNPAnnot.outs/Tcell/"
datasets <- read.table("../datasets.txt")$V1
for (i in 16:20){
   ##
   ii <- datasets[i]
   cat(ii, "\n")

   ## response motif 
   fn <- paste("/nfs/rprdata/julong/sc-atac/analyses.2021-02-05/3_motif/2_motif.activities.outs/5_4_",
      ii, ".motif.txt", sep="")
   motifList <- read.table(fn, header=F)$V1
   cat("motif", length(motifList), "\n") 
   snpmotif <- lapply(motifList, function(motif){
###
## cat(motif, "\n") 
      fn <- paste("./annot_jaspar2020/allsnp_", motif, ".bed.gz", sep="")   
      snp.motif <- try(fread(fn, header=F, data.table=F, stringsAsFactors=F), silent=T)
      if ( class(snp.motif)!="try-error"){ 
         snp.motif <- snp.motif%>%mutate(chr_pos=paste(V1, V2, sep=":"))
         xx <- snp.motif$chr_pos
      }else{
         xx <- NA
      }  
      xx  
   })
   snpmotif <- unlist(snpmotif[!is.na(snpmotif)])
   snpmotif <- unique(snpmotif)
   response.motif <- allsnp%>%filter(chr_pos%in%snpmotif)%>%dplyr::pull(V3)
   ###
   anno2 <- anno
   anno2[anno2[,2]!=0&anno2[,1]%in%response.motif,2] <- 3
    
   ## anno2 <- anno2%>%filter(SNP%in%testsnp) 
   ###
   opfn <- gzfile(paste(outdir, ii, "_torus.annot.gz", sep=""))
   write.table(anno2, opfn, sep="\t", row.names=F, col.names=T, quote=F)
   ###
}    

###
summ <- map_dfr(16:20, function(i){
    ##
    ii <- datasets[i]
    fn <-  paste(outdir, ii, "_torus.annot.gz", sep="")
    x <- fread(fn, header=T, data.table=F)
    dfsnp <- data.frame(condition=ii, "snp.1"=sum(x[,2]==1),"snp.2"= sum(x[,2]==2),"snp.3"= sum(x[,2]==3))
    dfsnp
})
###
opfn2 <- paste(outdir, "SNP_summary.csv", sep="")
write.csv(summ, opfn2, row.names=F)




#####################################
### response motif with direction ###
#####################################

rm(list=ls())

fn <- "SCAIP1-6_infor.txt.gz"
allsnp <- fread(fn, header=F, data.table=F)%>%mutate(chr_pos=paste(V1, V2, sep=":"))

###
testsnp <- read.table("../eQTL/eQTL_results/zzz_SNPlist.txt")$V1


###
### annotation
oneMCl <- "Tcell"
fn <-"./3_SNPAnnot.outs/Tcell/Tcell_Active_peak_torus.annot.gz"
anno <- fread(fn, header=T, data.table=F)
anno <- anno%>%filter(SNP%in%testsnp)


###
### cell-type specific motif 
fn <- "/nfs/rprdata/julong/sc-atac/analyses.2021-02-05/3_motif/1.3_motif.outs/2_4_Tcell.selectMotif.txt" 
motifList <- read.table(fn, header=T)[,1]
snpmotif <- lapply(motifList, function(motif){
###
## cat(motif, "\n") 
   fn <- paste("./annot_jaspar2020/allsnp_", motif, ".bed.gz", sep="")   
   snp.motif <- try(fread(fn, header=F, data.table=F, stringsAsFactors=F), silent=T)
   if ( class(snp.motif)!="try-error"){ 
      snp.motif <- snp.motif%>%mutate(chr_pos=paste(V1, V2, sep=":"))
      xx <- snp.motif$chr_pos
   }else{
      xx <- NA
   }  
   xx  
})
snpmotif <- unlist(snpmotif[!is.na(snpmotif)])
snpmotif <- unique(snpmotif)
snpmotif2 <- allsnp%>%filter(chr_pos%in%snpmotif)%>%dplyr::pull(V3)
snpmotif2 <- unique(snpmotif2)

anno[anno[,2]==1&anno[,1]%in%snpmotif2, 2] <- 2 


###
### response motifs
outdir <- "./3_SNPAnnot.outs/Tcell/"
datasets <- read.table("../datasets.txt")$V1
for (i in 16:20){
   ##
   ii <- datasets[i]
   cat(ii, "\n")

   ## response motif 
   fn <- paste("/nfs/rprdata/julong/sc-atac/analyses.2021-02-05/3_motif/2_motif.activities.outs/6_4_",
      ii, ".motif.txt", sep="")
   motifList <- read.table(fn, header=F)$V1
   cat("motif", length(motifList), "\n") 
   snpmotif <- lapply(motifList, function(motif){
###
## cat(motif, "\n") 
      fn <- paste("./annot_jaspar2020/allsnp_", motif, ".bed.gz", sep="")   
      snp.motif <- try(fread(fn, header=F, data.table=F, stringsAsFactors=F), silent=T)
      if ( class(snp.motif)!="try-error"){ 
         snp.motif <- snp.motif%>%mutate(chr_pos=paste(V1, V2, sep=":"))
         xx <- snp.motif$chr_pos
      }else{
         xx <- NA
      }  
      xx  
   })
   snpmotif <- unlist(snpmotif[!is.na(snpmotif)])
   snpmotif <- unique(snpmotif)
   response.motif <- allsnp%>%filter(chr_pos%in%snpmotif)%>%dplyr::pull(V3)
   ###
   anno2 <- anno
   anno2[anno2[,2]!=0&anno2[,1]%in%response.motif,2] <- 3
    
   ## anno2 <- anno2%>%filter(SNP%in%testsnp) 
   ###
   opfn <- gzfile(paste(outdir, "2_direction_", ii, "_torus.annot.gz", sep=""))
   write.table(anno2, opfn, sep="\t", row.names=F, col.names=T, quote=F)
   ###
}    

###
summ <- map_dfr(16:20, function(i){
    ##
    ii <- datasets[i]
    fn <-  paste(outdir,"2_direction_", ii, "_torus.annot.gz", sep="")
    x <- fread(fn, header=T, data.table=F)
    dfsnp <- data.frame(condition=ii, "snp.1"=sum(x[,2]==1),"snp.2"= sum(x[,2]==2),"snp.3"= sum(x[,2]==3))
    dfsnp
})
###
opfn2 <- paste(outdir, "SNP_summary.csv", sep="")
write.csv(summ, opfn2, row.names=F)    



###############################
### union of response motif ###
###############################

## response motif 
fn <- "/nfs/rprdata/julong/sc-atac/analyses.2021-02-05/3_motif/2_motif.activities.outs/7_4_Tcell_union.motif.txt"
###
motifList <- read.table(fn, header=F)$V1
cat("motif", length(motifList), "\n")

###
snpmotif <- lapply(motifList, function(motif){
###
## cat(motif, "\n") 
      fn <- paste("./annot_jaspar2020/allsnp_", motif, ".bed.gz", sep="")   
      snp.motif <- try(fread(fn, header=F, data.table=F, stringsAsFactors=F), silent=T)
      if ( class(snp.motif)!="try-error"){ 
         snp.motif <- snp.motif%>%mutate(chr_pos=paste(V1, V2, sep=":"))
         xx <- snp.motif$chr_pos
      }else{
         xx <- NA
      }  
      xx  
})
###
snpmotif <- unlist(snpmotif[!is.na(snpmotif)])
snpmotif <- unique(snpmotif)
response.motif <- allsnp%>%filter(chr_pos%in%snpmotif)%>%dplyr::pull(V3)

###
anno2 <- anno
anno2[anno2[,2]!=0&anno2[,1]%in%response.motif,2] <- 3
    
## anno2 <- anno2%>%filter(SNP%in%testsnp) 
###
opfn <- gzfile(paste(outdir, "3_union_Tcell_torus.annot.gz", sep=""))
write.table(anno2, opfn, sep="\t", row.names=F, col.names=T, quote=F)









