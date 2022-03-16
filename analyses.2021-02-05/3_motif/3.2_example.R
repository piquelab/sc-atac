##
library(Matrix)
library(tidyverse)
## library(parallel)
## library(data.table)
## library(future)
## library(purrr)
## library(furrr)
## library("BiocParallel")
## library(Rcpp)
## library(reshape)
library(qqman)
library(qvalue)
##
library(DESeq2)
library(biobroom)
library(ashr)
library(GenomicRanges)
library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(Signac)
library(clusterProfiler)
library(org.Hs.eg.db)
library(annotables)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(EnsDb.Hsapiens.v75)
library(ChIPseeker,lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
###
library(ggplot2)
library(cowplot,lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
library(grid)
library(gridExtra)
library(ggExtra)
library(RColorBrewer)
library(ggsci)
library(RColorBrewer)
library(viridis)
theme_set(theme_grey())

rm(list=ls())

outdir <- "./3.2_Example.outs/"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)


###
### read atac data
atac <- read_rds("../3_motif/1_motif.outs/1_scATAC.motif.rds")
## x <- atac@meta.data
## x$treat <- gsub(".*-ATAC-|_.*", "", x$NEW_BARCODE)
## atac <- AddMetaData(atac,x)

###
### peak annotation
fn <- "../2_Differential/2.2_compareRNAandATAC.outs/2_annot.ChIPseeker.rds"
anno <- read_rds(fn)%>%as.data.frame()%>%
       mutate(peak=paste(gsub("chr","",seqnames), start, end, sep="-"))


####
#### fun-1, get gene id
getGene <- function(atac, anno, motif.name, resDP){

   ### motif data   
   motif <- Motifs(atac)
   x <- motif@data
   motif.ID <- ConvertMotifID(motif, motif.name)

   ### peak containing motifs 
   peaks <- rownames(x)
   peak.sel <- peaks[x[,motif.ID[1]]==1] ## peak containing motifs

  ###
   peak.final <- intersect(peak.sel, unique(resDP$gene))  
   anno2 <- anno%>%dplyr::filter(peak%in%peak.final, abs(distanceToTSS)<100000) ##100000

   ### gene list in a set of peaks that contain specific motif and are differentially accessible
   gene <- unique(anno2$geneId)
   gene 
}    


###
### fun-2, average gene expression
avePathway <- function(X){

   ### filtering more missing value
   ii <- apply(!is.na(X), 1, sum)
   X <- X[ii>20,]
   bti <- colnames(X)
   cvt <- str_split(bti, "_", simplify=T)
   cvt <- data.frame(rn=bti, MCls=cvt[,1],treats=cvt[,2],sampleID=cvt[,3],Batch=cvt[,4])%>%
       mutate(comb=paste(MCls, Batch, sep="_"))
   comb <- unique(cvt$comb)
   for (ii in comb){
      bti0 <- cvt%>%dplyr::filter(comb==ii)%>%dplyr::pull(rn)
      x <- X[,bti0]
      x.mean <- apply(x, 1, mean, na.rm=T)
      x.scale <- sweep(x, 1, x.mean, "-")
      X[,bti0] <- x.scale
   }
   pathway <- apply(X, 2, mean, na.rm=T)
}

###
### fun-3
getData <- function(gene, datatype="bulk"){

   ### bulk 
   if ( datatype=="bulk"){
      fn <- "/nfs/rprdata/julong/SCAIP/analyses/SCAIP-B1-6_2020.03.23/6_DEG.CelltypeNew_output/Filter2/YtX_sel.comb.RData"
      load(fn)
      rn <- gsub("\\..*", "", rownames(YtX_sel))
      rownames(YtX_sel) <- rn
      X <- YtX_sel[rn%in%gene,]+1
      counts <- colSums(YtX_sel)
      X <- sweep(X, 2, counts, "/")
      X <- X*1e+06
      X <- log2(X)
     
      bti <- colnames(X)
      cvt <- str_split(bti, "_", simplify=T)
      cvt <- data.frame(rn=bti, MCls=cvt[,1], treats=gsub("-EtOH", "", cvt[,2]),
                        sampleID=cvt[,3], Batch=cvt[,4])
      cvt$y <- avePathway(X)
   ##
   } ### end bulk
    
    
   #### variability    
   if ( datatype=="NB.phi"){
      fn <- "/nfs/rprdata/julong/SCAIP/analyses/SCAIP-B1-6_2020.03.23/10_RNA.Variance_output/tmp10/1.2_Sel.PhxNew.RData"
      load(fn)
      rn <- gsub("\\..*", "", rownames(PhxNew2))
      rownames(PhxNew2) <- rn

      X <- PhxNew2[rn%in%gene,]
      X <- log2(X)
      bti <- colnames(PhxNew2)
      cvt <- str_split(bti, "_", simplify=T)
      cvt <- data.frame(rn=bti, MCls=cvt[,1], treats=gsub("-EtOH", "", cvt[,2]), sampleID=cvt[,3], Batch=cvt[,4])
      cvt$y <- avePathway(X)
   } ### end variability

    
   MCls <- sort(unique(cvt$MCls))
   ### 
   cvt <- map_dfr(MCls, function(oneMCl){
      cvt2 <- cvt%>%dplyr::filter(MCls==oneMCl)
      df0 <- cvt2%>%dplyr::filter(treats=="CTRL")%>%
         dplyr::select(sampleID, y)%>%
         dplyr::rename("y0"="y")
      cvt2 <- cvt2%>%left_join(df0, by=c("sampleID"))%>%mutate(y2=y-y0)
      cvt2
   })
   ###    
   cvt
}
    


## lab1 <- c("CTRL"="CTRL",
lab1 <- c("LPS"="LPS", "LPS-DEX"="LPS+DEX",
   "PHA"="PHA", "PHA-DEX"="PHA+DEX")
## col1 <- c("CTRL"="#828282",
col1 <- c("LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
   "PHA"="#a6cee3", "PHA-DEX"="#1f78b4")


## motif.name <- "RELA"
res <- read_rds("../2_Differential/1.2_DiffPeak.outs/2.0_DESeq.results.rds")%>%
    as.data.frame()

motifList <- c("NR3C1", "NR3C2", "RELA", "FOS::JUN", "NFKB1", "NFKB2")


##################
### expression ###
##################

## NR3C2
###
###
for (i in 1:2){
   ##
   motif <- motifList[i]
   cat(i, motif, "\n") 
 ### differential peaks 
   res2 <- res%>%
      dplyr::filter(p.adjusted<0.1, abs(estimate)>0.5, estimate>0, MCls=="Tcell", grepl("DEX",contrast))
    
   gene <- getGene(atac, anno, motif.name=motif, resDP=res2)
   cvt <- getData(gene=gene, datatype="bulk")
###
###
   cvt2 <- cvt%>%dplyr::filter(treats!="CTRL", MCls=="Tcell", Batch=="SCAIP6")

   p <- ggplot(cvt2,aes(x=factor(treats), y=y, fill=treats))+
      geom_boxplot(outlier.size=0.8)+
      ylab("Expression value")+
      scale_y_continuous(expand=expansion(mult=c(0, 0.2)))+
      scale_fill_manual("", values=col1, labels=lab1)+
      scale_x_discrete("", labels=lab1)+
      ggtitle(bquote(~italic(.(motif))~"(T cell)"))+
      theme_bw()+
      theme(axis.text.x=element_text(angle=-90, hjust=0, size=10),
          axis.title.x=element_blank(),
          axis.title.y=element_text(size=12),
          plot.title=element_text(hjust=0.5, size=12),
          legend.position="none")

###
   figfn <- paste("./3.2_Example.outs/Figure", i, "_", motif, ".boxplot.png", sep="")
   png(figfn, width=480, height=420, res=120)
   print(p)
   dev.off()
###
}##

###
###
for (i in 3:5){
   motif <- motifList[i]
   cat(i, motif, "\n")
 ### differential peaks 
   res2 <- res%>%
      dplyr::filter(p.adjusted<0.1, abs(estimate)>0.5, estimate<0, MCls=="Tcell", grepl("DEX",contrast))
    
   gene <- getGene(atac, anno, motif.name=motif, resDP=res2)
   cvt <- getData(gene=gene, datatype="bulk")
###
###
   cvt2 <- cvt%>%dplyr::filter(treats!="CTRL", MCls=="Tcell", Batch=="SCAIP6")

   p <- ggplot(cvt2,aes(x=factor(treats), y=y, fill=treats))+
      geom_boxplot(outlier.size=0.8)+
      ylab("Expression value")+
      scale_y_continuous(expand=expansion(mult=c(0, 0.2)))+
      scale_fill_manual("", values=col1, labels=lab1)+
      scale_x_discrete("", labels=lab1)+
      ggtitle(bquote(~italic(.(motif))~"(T cell)"))+
      theme_bw()+
      theme(axis.text.x=element_text(angle=-90, hjust=0, size=10),
          axis.title.x=element_blank(),
          axis.title.y=element_text(size=12),
          plot.title=element_text(hjust=0.5, size=12),
          legend.position="none")

###
   figfn <- paste("./3.2_Example.outs/Figure", i, "_", motif, ".boxplot.png", sep="")
   png(figfn, width=480, height=420, res=120)
   print(p)
   dev.off()
}



###################
### variability ###
###################

for (i in 1:2){
   ##
   motif <- motifList[i]
   cat(i, motif, "\n") 
 ### differential peaks 
   res2 <- res%>%
      dplyr::filter(p.adjusted<0.1, abs(estimate)>0.5, estimate>0, MCls=="Tcell", grepl("DEX",contrast))
    
   gene <- getGene(atac, anno, motif.name=motif, resDP=res2)
   cvt <- getData(gene=gene, datatype="NB.phi")
###
###
   cvt2 <- cvt%>%dplyr::filter(treats!="CTRL", MCls=="Tcell", Batch=="SCAIP6")

   p <- ggplot(cvt2,aes(x=factor(treats), y=y, fill=treats))+
      geom_boxplot(outlier.size=0.8)+
      ylab("Variability")+
      scale_y_continuous(expand=expansion(mult=c(0, 0.2)))+
      scale_fill_manual("", values=col1, labels=lab1)+
      scale_x_discrete("", labels=lab1)+
      ggtitle(bquote(~italic(.(motif))~"(T cell)"))+
      theme_bw()+
      theme(axis.text.x=element_text(angle=-90, hjust=0, size=10),
          axis.title.x=element_blank(),
          axis.title.y=element_text(size=12),
          plot.title=element_text(hjust=0.5, size=12),
          legend.position="none")

###
   figfn <- paste("./3.2_Example.outs/Figure", i, "_", motif, "_va.boxplot.png", sep="")
   png(figfn, width=480, height=420, res=120)
   print(p)
   dev.off()
###
}##




### nearby gene expression
## gene0 <- res5$geneId[1]
## symbol0 <- res5$SYMBOL[1]
## cvt <- getData(gene=gene0, datatype="RNA")
## ###
## ###
## cvt2 <- cvt%>%dplyr::filter(treats!="CTRL", MCls=="Tcell")
## p <- ggplot(cvt2,aes(x=factor(treats), y=yscale2, fill=treats))+
##    geom_boxplot(outlier.size=0.8)+
##    ylab(bquote(~log[2]~"(Expression)"))+
##    scale_y_continuous(expand=expansion(mult=c(0, 0.2)))+
##    scale_fill_manual("", values=col1, labels=lab1)+
##    scale_x_discrete("", labels=lab1)+
##    ggtitle(bquote(~italic(.(symbol0))~"(T cell)"))+
##    theme_bw()+
##    theme(axis.text.x=element_text(angle=-90, hjust=0, size=10),
##          axis.title.x=element_blank(),
##          axis.title.y=element_text(size=12),
##          plot.title=element_text(hjust=0.5, size=12),
##          legend.position="none")

## ###
## figfn <- paste("./5_Example.outs/", motif.name, "/Figure", i, ".3_", symbol0, ".boxplot.png", sep="")
## png(figfn, width=480, height=420, res=120)
## print(p)
## dev.off()

## ### gene variability
## cvt <- getData(gene=gene0, datatype="NB.phi")
## cvt2 <- cvt%>%dplyr::filter(treats!="CTRL", MCls=="Tcell")
## p <- ggplot(cvt2,aes(x=factor(treats), y=yscale2, fill=treats))+
##    geom_boxplot(outlier.size=0.8)+
##    ylab(bquote(~log[2]~"(Variability)"))+
##    scale_y_continuous(expand=expansion(mult=c(0, 0.2)))+
##    scale_fill_manual("", values=col1, labels=lab1)+
##    scale_x_discrete("", labels=lab1)+
##    ggtitle(bquote(~italic(.(symbol0))~"(T cell)"))+
##    theme_bw()+
##    theme(axis.text.x=element_text(angle=-90, hjust=0, size=10),
##          axis.title.x=element_blank(),
##          axis.title.y=element_text(size=12),
##          plot.title=element_text(hjust=0.5, size=12),
##          legend.position="none")

## ###
## figfn <- paste("./5_Example.outs/", motif.name, "/Figure", i, ".4_", symbol0, ".va.boxplot.png", sep="")
## png(figfn, width=480, height=420, res=120)
## print(p)
## dev.off()



## geneID <- bitr(geneID=unique(resDE$gene), fromType="ENSEMBL", toType="SYMBOL",
##    OrgDb=org.Hs.eg.db)
## geneID <- geneID%>%dplyr::rename("gene"="ENSEMBL")

## resDE <- resDE%>%dplyr::left_join(geneID, by="gene")
