###
###
library(tidyverse)
library(data.table)
## library(Seurat)
## library(Signac)
library(Matrix)
 
### annotation SNP hit by peak, cell-type motif and response motif
### grid of pct values
### edit by Julong wei, 2-14-2022


### 0
### check motif if exists
## motifList <- read.table("motifList2020.txt")$V1
## for (i in 1:length(motifList)){
## ##
##    motif <- motifList[i]  
##    fn <- paste("./annot_jaspar2020/allsnp_", motif, ".bed.gz", sep="")   
##    snp.motif <- try(fread(fn, header=F, data.table=F, stringsAsFactors=F), silent=T)
##    ##
    
##    ## if ( class(snp.motif)=="try-error"){
##    ##    cat(i, motif, "\n")
##    ##  }
##    if ( !file.exists(fn)|file.size(fn)==0){
##        cat(i, motif, "\n")
##     }   
    
## }


### 1
### generate cell-type active peak bed file

prefix.in <- "/nfs/rprdata/julong/sc-atac/analyses.2021-02-05/3_motif/1.3_motif.outs/"

 
pct_grid <- c(0.05, 0.02, 0.01)
MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
###
for (i in 1:4){
   ##
   pct0 <- pct_grid[2]   
   outdir <- paste("./4_SNPAnnot.outs/pct_", pct0, "/",  sep="")
   if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)
    
   ###
   fn <- paste(prefix.in, "pct_", pct0, "/1_cell-type_active.peaks.rds", sep="")
   peak <- read_rds(fn)
   ## 
   ii <- i+1
   peak2 <- peak[peak[,ii]==1,] 
   peakSel <- peak2$peakAll
    
   bed <- str_split(peakSel, "-", simplify=T)
    
   ###
   opfn <- paste(outdir, MCls[i], "_Active_peak.bed", sep="")
   write.table(bed, opfn, sep="\t", row.names=F, quote=F, col.names=F)
   ##
   cat(MCls[i], nrow(bed), "\n")  
}


### 2
### test if SNP is hit by cell-type active peak
## 4.2_intersect.sh


### 3
### torus format
rm(list=ls())

fn <- "SCAIP1-6_infor.txt.gz"
allsnp <- fread(fn, header=F, data.table=F)%>%mutate(chr_pos=paste(V1, V2, sep=":"))

MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
pct_grid <- c(0.05, 0.02, 0.01)
## for(i in 1:length(pct_grid)){
pct0 <- pct_grid[2]
###
for (i in 1:4){
   ## 
   outdir <- paste("4_SNPAnnot.outs/pct_", pct0, "/", sep="")
   fn <- paste(outdir, MCls[i], "_Active_peak_snp.bed.gz", sep="")
   snpPeak <- fread(fn, header=F, data.table=F)
   snpPeak2 <- snpPeak%>%mutate(chr_pos=paste(V1, V2, sep=":"))%>%dplyr::filter(V4>0)%>%dplyr::pull(chr_pos)
   ## 
   peaking_d <- ifelse(allsnp$chr_pos%in%snpPeak2, 1, 0)
   anno <- data.frame("SNP"=allsnp$V3, "peaking_d"=peaking_d)
   ##
   opfn <- gzfile(paste(outdir,  MCls[i], "_Active_peak_torus.annot.gz", sep=""))
   write.table(anno, opfn, sep="\t", row.names=F, col.names=T, quote=F) 
   ##
   cat(MCls[i], sum(anno[,2]), "\n") 
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

#############################################################################
### 4, add SNP annotation hit by peak, T-cell motifs and treatment motifs ###
#############################################################################

### definition of response motifs in the (1) way, for treatment respectively, not consider direction

rm(list=ls())


fn <- "SCAIP1-6_infor.txt.gz"
allsnp <- fread(fn, header=F, data.table=F)%>%mutate(chr_pos=paste(V1, V2, sep=":"))

###
testsnp <- read.table("../eQTL/eQTL_results/zzz_SNPlist.txt")$V1


###
### annotation
pct_grid <- c(0.05, 0.02, 0.01)
for ( k in 1:length(pct_grid)){
###
   pct0 <- pct_grid[k]
   cat("pct", pct0)
   outdir <- paste("./4_SNPAnnot.outs/pct_", pct0, "/", sep="") 
   indir <- paste("/nfs/rprdata/julong/sc-atac/analyses.2021-02-05/3_motif/1.3_motif.outs/pct_", pct0, "/", sep="") 
   ##    
   fn <- paste(outdir,  "Tcell_Active_peak_torus.annot.gz", sep="")    
   anno <- fread(fn, header=T, data.table=F)
   anno <- anno%>%filter(SNP%in%testsnp)

###
### cell-type specific motif 
   fn <- paste(indir, "2_4_Tcell.selectMotif.txt", sep="") 
   motifList <- read.table(fn, header=T)[,1]
   cat("cell-type motif", length(motifList), "\n")
   ##    
   snpmotif <- lapply(motifList, function(motif){
###
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
snpmotif2 <- allsnp%>%filter(chr_pos%in%snpmotif)%>%dplyr::pull(V3)
snpmotif2 <- unique(snpmotif2)

anno[anno[,2]==1&anno[,1]%in%snpmotif2, 2] <- 2 


###
### response motifs
datasets <- read.table("../datasets.txt")$V1
for (i in 16:20){
   ##
   ii <- datasets[i]
   ## cat(ii, "\n")

   ## response motif 
   fn <- paste("/nfs/rprdata/julong/sc-atac/analyses.2021-02-05/3_motif/2_motif.activities.outs/5_4_",
      ii, ".motif.txt", sep="")
   motifList <- read.table(fn, header=F)$V1
   ## cat("motif", length(motifList), "\n") 
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
} ### End conditions

    
} ## End

    
###
summ2 <- map_dfr(1:length(pct_grid), function(k){
    pct0 <- pct_grid[k]
###
summ <- map_dfr(16:20, function(i){
    ##
    ii <- datasets[i]
    outdir <- paste("./4_SNPAnnot.outs/pct_", pct0, "/", sep="")
    fn <-  paste(outdir, ii, "_torus.annot.gz", sep="")
    x <- fread(fn, header=T, data.table=F)
    dfsnp <- data.frame(condition=ii, "snp.1"=sum(x[,2]==1),
       "snp.2"= sum(x[,2]==2), "snp.3"= sum(x[,2]==3), "snp.total"=sum(x[,2]!=0), "pct"=pct0)
    dfsnp
})
summ
})

###
opfn2 <- "./4_SNPAnnot.outs/SNP_summary.csv"
write.csv(summ2, opfn2, row.names=F)



#############################################################################
### 5, add SNP annotation hit by peak, cell-type motifs and treatment motifs ###
#############################################################################

### definition of response motifs in the (3) way, union of response motifs

###############################
### union of response motif ###
###############################

rm(list=ls())

###
###

fn <- "SCAIP1-6_infor.txt.gz"
allsnp <- fread(fn, header=F, data.table=F)%>%mutate(chr_pos=paste(V1, V2, sep=":"))
###
## testsnp <- read.table("../eQTL/eQTL_results/zzz_SNPlist.txt")$V1

##
###
###

pct_grid <- c(0.05, 0.02, 0.01)
MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
pct0 <- pct_grid[2]

###
### for each cell-type 
for ( i in 1:4){

   oneMCl <- MCls[i]    
  
   ###
   ### SNP hit by response motif 
   fn <- paste("/nfs/rprdata/julong/sc-atac/analyses.2021-02-05/3_motif/2_motif.activities.outs/7.", i,
      "_", oneMCl, "_union.motif.txt", sep="")         
   ### 
   motifList <- read.table(fn, header=T)$gene
   cat("motif", oneMCl, length(motifList), "\n")

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
   ### annotation
   outdir <- paste("./4_SNPAnnot.outs/pct_", pct0, "/", sep="")
   ##    
   fn <- paste(outdir,  oneMCl, "_Active_peak_torus.annot.gz", sep="")    
   anno <- fread(fn, header=T, data.table=F)
   ## anno <- anno%>%filter(SNP%in%testsnp)

   ###
   ### cell-type specific motif 
   fn <- paste("/nfs/rprdata/julong/sc-atac/analyses.2021-02-05/3_motif/1.3_motif.outs/pct_", pct0, "/2.",
      i, "_", oneMCl, ".selectMotif.txt", sep="") 
   motifList <- read.table(fn, header=T)[,1]
   cat("cell-type motif", oneMCl, length(motifList), "\n")
   ##    
   snpmotif <- lapply(motifList, function(motif){
###
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
   snpmotif2 <- allsnp%>%filter(chr_pos%in%snpmotif)%>%dplyr::pull(V3)
   snpmotif2 <- unique(snpmotif2)
   ## cell-type motifs
   anno[anno[,2]==1&anno[,1]%in%snpmotif2, 2] <- 2 

   ### response motif
   anno2 <- anno
   anno2[anno2[,2]!=0&anno2[,1]%in%response.motif,2] <- 3
   ### 
   cat(table(anno2[,2]), "\n")
   ### output    
   opfn <- gzfile(paste(outdir, "3_", oneMCl, "_union_torus.annot.gz", sep=""))
   write.table(anno2, opfn, sep="\t", row.names=F, col.names=T, quote=F)
}
### End


## ###
## pct_grid <- c(0.05, 0.02, 0.01)
## summ2 <- map_dfr(1:length(pct_grid), function(k){
##     pct0 <- pct_grid[k]
##     ##
##     outdir <- paste("./4_SNPAnnot.outs/pct_", pct0, "/", sep="")
##     fn <-  paste(outdir, "3_union_Tcell_torus.annot.gz", sep="")
##     x <- fread(fn, header=T, data.table=F)
##     dfsnp <- data.frame("snp.1"=sum(x[,2]==1),
##        "snp.2"= sum(x[,2]==2), "snp.3"= sum(x[,2]==3), "snp.total"=sum(x[,2]!=0), "pct"=pct0)
##     dfsnp
## })
## ###

## pct0 <- 0.02





