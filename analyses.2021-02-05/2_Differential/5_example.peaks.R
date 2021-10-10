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

outdir <- "./5_Example.outs/"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)


###
### read atac data
atac <- read_rds("../3_motif/1_motif.outs/1_scATAC.motif.rds")
x <- atac@meta.data
x$treat <- gsub(".*-ATAC-|_.*", "", x$NEW_BARCODE)
atac <- AddMetaData(atac,x)

###
### pick peaks containing motifs
motif <- Motifs(atac)
x <- motif@data

motif.name <- c("NR3C1","NR3C2")
## motif.name <- "RELA"
motif.ID <- ConvertMotifID(motif, motif.name)

## dir.create("./5_Example.outs/RELA/", showWarnings=F, recursive=F)

peaks <- rownames(x)
peak.sel <- peaks[x[,motif.ID[1]]==1] ## peak containing motifs


################################################
### Differential expressed or variable genes ###
################################################

fn <- "/nfs/rprdata/julong/SCAIP/analyses/SCAIP-B1-6_2020.03.23/6_DEG.CelltypeNew_output/Filter2/2_meta.rds"
resDE <- read_rds(fn)%>%as.data.frame()
resDE2 <- resDE%>%dplyr::filter(qval<0.1,abs(beta)>0.5, MCls=="Tcell") 
DEG <- unique(resDE2$rn)

fn <- "/nfs/rprdata/julong/SCAIP/analyses/SCAIP-B1-6_2020.03.23/10_RNA.Variance_output/tmp9/3_phiNew.meta"
resDV <- read.table(fn, header=T)%>%as.data.frame()
resDV2 <- resDV%>%dplyr::filter(qval<0.1,abs(beta)>0.5, MCls=="Tcell") 
DVG <- unique(resDV2$rn)




#######################################################
### Differential peaks that contain specific motif  ###
### and annotate peaks using closest genes          ###
#######################################################

res <- read_rds("./1.2_DiffPeak.outs/2.0_DESeq.results.rds")%>%
    as.data.frame()%>%
   dplyr::filter(p.adjusted<0.1, abs(estimate)>0.5)

res2 <- res%>%dplyr::filter(gene%in%peak.sel)


### annotation
anno <- read_rds("./2.2_compareRNAandATAC.outs/2_annot.ChIPseeker.rds")%>%
   as.data.frame%>%
   mutate(peak=paste(gsub("chr","",seqnames), start, end, sep="-"))

anno2 <- anno%>%dplyr::select(seqnames, start, end, geneId, SYMBOL, peak)

### 
res3 <- res2%>%left_join(anno2, by=c("gene"="peak"))%>%
   arrange(desc(abs(statistic)))

res4 <- res3%>%dplyr::filter(MCls=="Tcell")

res4 <- res4%>%mutate(rn=paste(MCls, contrast, geneId, sep="_"),
   is_DEG=ifelse(rn%in%DEG, 1, 0), is_DVG=ifelse(rn%in%DVG,1,0))

res5 <- res4%>%dplyr::filter(is_DEG==1,estimate<0)
###
### plot specific region that was differentially expressed and contain motif and DEG
atac2 <- subset(atac,subset=MCls=="Tcell")

i <- 1
peak0 <- res5$gene[1]
ranges.show <- StringToGRanges(peak0)
p <- CoveragePlot(atac2,
   region=peak0, region.highlight=ranges.show,
   group.by="treat",
   extend.upstream=8e+03, extend.downstream=8e+03)&
   scale_fill_manual(values=c("CTRL"="#828282",
      "LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
      "PHA"="#a6cee3", "PHA-DEX"="#1f78b4"))    

figfn <- paste("./5_Example.outs/", motif.name, "/Figure", i, ".1_", peak0,
   ".coverage.png", sep="")
png(figfn, width=580, height=400, res=120)
print(p)
dev.off()

###
### mancol


###############
### boxplot ###
###############

###
adjGene <- function(cvt, center=T, contrast=T){
    
   ###centralization 
   cvt <- cvt%>%mutate(comb=paste(MCls, Batch, sep="_"))
   if(center){
      cvt <- cvt%>%group_by(comb)%>%mutate(yscale=y-mean(y,na.rm=T))%>%ungroup()
   }else{
      cvt <- cvt%>%mutate(yscale=y)
   }
    
   ### minus contrl group
   if (contrast){ 
      cvt0 <- cvt%>%dplyr::filter(treats=="CTRL")
      y0 <- cvt0$yscale
      names(y0) <- cvt0$sampleID
      y0[is.na(y0)] <- 0
      cvt$y0 <- y0[cvt$sampleID]
      cvt <- cvt%>%mutate(yscale2=yscale-y0)
   }
    
   cvt    
}


###
getData <- function(gene, datatype="ATAC"){

   if ( datatype=="ATAC"){  
      fn <- "./1.2_DiffPeak.outs/1_YtX.sel.rds"
      YtX <- read_rds(fn)
###
      bti <- colnames(YtX)
      cvt <- str_split(bti, "_", simplify=T)
      cvt <- data.frame(rn=bti, MCls=cvt[,1], treats=cvt[,2],
                        sampleID=cvt[,3], Batch="SCAIP6")   
##
      rn <- rownames(YtX)
      X <- YtX[rn%in%gene,]+1   
      counts <- colSums(YtX)

      X <- (X/counts)*1e+06
      cvt$y <- log2(X)
   }
                                                            
    
   if ( datatype=="RNA"){
      fn <- "/nfs/rprdata/julong/SCAIP/analyses/SCAIP-B1-6_2020.03.23/6_DEG.CelltypeNew_output/Filter2/YtX_sel.comb.RData"
      load(fn)
      rn <- gsub("\\..*", "", rownames(YtX_sel))
      rownames(YtX_sel) <- rn
      X <- YtX_sel[rn%in%gene,]+1
      counts <- colSums(YtX_sel)
      X <- ( X/counts)*1e+06

      bti <- colnames(YtX_sel)
      cvt <- str_split(bti, "_", simplify=T)
      cvt <- data.frame(rn=bti, MCls=cvt[,1], treats=gsub("-EtOH", "", cvt[,2]),
                        sampleID=cvt[,3], Batch=cvt[,4])
      cvt$y <- log2(X)
  }

  if ( datatype=="NB.phi"){
     fn <- "/nfs/rprdata/julong/SCAIP/analyses/SCAIP-B1-6_2020.03.23/10_RNA.Variance_output/tmp10/1.2_Sel.PhxNew.RData"
     load(fn)
     rn <- gsub("\\..*", "", rownames(PhxNew2))
     rownames(PhxNew2) <- rn

     X <- PhxNew2[rn%in%gene,]
     bti <- colnames(PhxNew2)
     cvt <- str_split(bti, "_", simplify=T)
     cvt <- data.frame(rn=bti, MCls=cvt[,1], treats=gsub("-EtOH", "", cvt[,2]), sampleID=cvt[,3], Batch=cvt[,4])
     cvt$y <- log2(X)
  }

  cvt <- adjGene(cvt)  
    
}


## lab1 <- c("CTRL"="CTRL",
lab1 <- c("LPS"="LPS", "LPS-DEX"="LPS+DEX",
   "PHA"="PHA", "PHA-DEX"="PHA+DEX")
## col1 <- c("CTRL"="#828282",
col1 <- c("LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
   "PHA"="#a6cee3", "PHA-DEX"="#1f78b4")

peak0 <- res5$gene[1]
cvt <- getData(gene=peak0, datatype="ATAC")
###
###
cvt2 <- cvt%>%dplyr::filter(treats!="CTRL", MCls=="Tcell")
p <- ggplot(cvt2,aes(x=factor(treats), y=yscale2, fill=treats))+
   geom_boxplot(outlier.size=0.8)+
   ylab(bquote(~log[2]~"(Expression)"))+
   scale_y_continuous(expand=expansion(mult=c(0, 0.2)))+
   scale_fill_manual("", values=col1, labels=lab1)+
   scale_x_discrete("", labels=lab1)+
   ggtitle(bquote(~.(peak0)~"(T cell)"))+
   theme_bw()+
   theme(axis.text.x=element_text(angle=-90, hjust=0, size=10),
         axis.title.x=element_blank(),
         axis.title.y=element_text(size=12),
         plot.title=element_text(hjust=0.5, size=12),
         legend.position="none")

###
figfn <- paste("./5_Example.outs/", motif.name, "/Figure", i, ".2_", peak0, ".boxplot.png", sep="")
png(figfn, width=480, height=420, res=120)
print(p)
dev.off()


### nearby gene expression
gene0 <- res5$geneId[1]
symbol0 <- res5$SYMBOL[1]
cvt <- getData(gene=gene0, datatype="RNA")
###
###
cvt2 <- cvt%>%dplyr::filter(treats!="CTRL", MCls=="Tcell")
p <- ggplot(cvt2,aes(x=factor(treats), y=yscale2, fill=treats))+
   geom_boxplot(outlier.size=0.8)+
   ylab(bquote(~log[2]~"(Expression)"))+
   scale_y_continuous(expand=expansion(mult=c(0, 0.2)))+
   scale_fill_manual("", values=col1, labels=lab1)+
   scale_x_discrete("", labels=lab1)+
   ggtitle(bquote(~italic(.(symbol0))~"(T cell)"))+
   theme_bw()+
   theme(axis.text.x=element_text(angle=-90, hjust=0, size=10),
         axis.title.x=element_blank(),
         axis.title.y=element_text(size=12),
         plot.title=element_text(hjust=0.5, size=12),
         legend.position="none")

###
figfn <- paste("./5_Example.outs/", motif.name, "/Figure", i, ".3_", symbol0, ".boxplot.png", sep="")
png(figfn, width=480, height=420, res=120)
print(p)
dev.off()

### gene variability
cvt <- getData(gene=gene0, datatype="NB.phi")
cvt2 <- cvt%>%dplyr::filter(treats!="CTRL", MCls=="Tcell")
p <- ggplot(cvt2,aes(x=factor(treats), y=yscale2, fill=treats))+
   geom_boxplot(outlier.size=0.8)+
   ylab(bquote(~log[2]~"(Variability)"))+
   scale_y_continuous(expand=expansion(mult=c(0, 0.2)))+
   scale_fill_manual("", values=col1, labels=lab1)+
   scale_x_discrete("", labels=lab1)+
   ggtitle(bquote(~italic(.(symbol0))~"(T cell)"))+
   theme_bw()+
   theme(axis.text.x=element_text(angle=-90, hjust=0, size=10),
         axis.title.x=element_blank(),
         axis.title.y=element_text(size=12),
         plot.title=element_text(hjust=0.5, size=12),
         legend.position="none")

###
figfn <- paste("./5_Example.outs/", motif.name, "/Figure", i, ".4_", symbol0, ".va.boxplot.png", sep="")
png(figfn, width=480, height=420, res=120)
print(p)
dev.off()



## geneID <- bitr(geneID=unique(resDE$gene), fromType="ENSEMBL", toType="SYMBOL",
##    OrgDb=org.Hs.eg.db)
## geneID <- geneID%>%dplyr::rename("gene"="ENSEMBL")

## resDE <- resDE%>%dplyr::left_join(geneID, by="gene")
