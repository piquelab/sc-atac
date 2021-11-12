###
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


####
outdir <- "./3.1_Example.outs/"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)


###
## getData <- function(gene, celltype){
##    fn <- "./1.2_DiffPeak.outs/1_YtX.sel.rds"
##    YtX <- read_rds(fn)
## ##
##    bti <- colnames(YtX)
##    cvt <- str_split(bti, "_", simplify=T)
##    cvt <- data.frame(rn=bti, MCls=cvt[,1], treats=cvt[,2], sampleID=cvt[,3])
   
## ##
##    rn <- rownames(YtX)
##    X <- YtX[rn%in%gene,]+1   
##    counts <- colSums(YtX)

##    X <- (X/counts)*1e+06
##    cvt$y <- log2(X)
##    cvt2 <- cvt%>%dplyr::filter(MCls==celltype)
##    cvt2
## }

getData <- function(motif){
   ###
   X <- read_rds("./2_motif.activities.outs/2_motif.ave.rds")
   bti <- colnames(X)
   cvt <- str_split(bti, "_", simplify=T)
   cvt <- data.frame(rn=bti, MCls=cvt[,1], treats=cvt[,2], sampleID=cvt[,3])
   cvt$y <- X[motif,]
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


##
getData2 <- function(atac, motif.name.sel){
###    
   X <- atac@assays$chromvar@data
   motif <- Motifs(atac)
   motif.name <- ConvertMotifID(object=motif, id=rownames(X))
   rownames(X) <- motif.name

   meta <- atac@meta.data%>%mutate(treats=gsub(".*-ATAC-|_.*", "", NEW_BARCODE))

   cvt <- data.frame(y=X[motif.name.sel,],treats=meta$treats, MCls=meta$MCls)
   cvt
}



lab1 <- c("CTRL"="CTRL",
   "LPS"="LPS", "LPS-DEX"="LPS+DEX",
   "PHA"="PHA", "PHA-DEX"="PHA+DEX")
col1 <- c("CTRL"="#828282",
   "LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
   "PHA"="#a6cee3", "PHA-DEX"="#1f78b4")

atac <- read_rds("./2_motif.activities.outs/1_scATAC.motifActivities.rds")


motifList <- c("NR3C1", "RELA", "NR3C2")
###
for ( i in 1:3){
    
motif <- motifList[i]
### average motif activities for each individual
cvt <- getData(motif)
###
cvt2 <- cvt%>%dplyr::filter(treats!="CTRL")    
p <- ggplot(cvt2, aes(x=factor(treats), y=y2, fill=treats, color=treats))+
   geom_boxplot(outlier.shape=NA)+
    geom_jitter(width=0.15, size=0.3)+ 
   ylab("motif activities")+
   scale_y_continuous(expand=expansion(mult=c(0, 0.2)))+
   scale_fill_manual("", values=col1, labels=lab1)+
   scale_colour_manual("",values=col1)+ 
   scale_x_discrete("", labels=lab1)+
   facet_wrap(~MCls, nrow=2, scales="free")+ 
   ggtitle(bquote(~italic(.(motif))))+
   theme_bw()+
   theme(axis.text.x=element_text(angle=-90, hjust=0, size=8),
         axis.title.x=element_blank(),
         axis.title.y=element_text(size=12),
         plot.title=element_text(hjust=0.5, size=12),
         strip.text=element_text(size=12),
         legend.position="none")

###
figfn <- paste("./3.1_Example.outs/Figure", i, ".1_", motif, ".boxplot.png", sep="")
png(figfn, width=500, height=500, res=120)
print(p)
dev.off()
}
## motif activities for each cell
## cvt2 <- getData2(atac, motif)%>%dplyr::filter(MCls!="DC")
## p <- ggplot(cvt2,aes(x=factor(treats), y=y, fill=treats))+
##    geom_boxplot(outlier.size=0.5)+
##    ylab("motif activities")+
##    scale_y_continuous(expand=expansion(mult=c(0, 0.2)))+
##    scale_fill_manual("", values=col1, labels=lab1)+
##    ## scale_colour_manual("", values=col1)+ 
##    scale_x_discrete("", labels=lab1)+
##    facet_wrap(~MCls, nrow=2, scales="free")+ 
##    ggtitle(bquote(~italic(.(motif))))+
##    theme_bw()+
##    theme(axis.text.x=element_text(angle=-90, hjust=0, size=8),
##          axis.title.x=element_blank(),
##          axis.title.y=element_text(size=12),
##          plot.title=element_text(hjust=0.5, size=12),
##          strip.text=element_text(size=12),
##          legend.position="none")

## ###
## figfn <- paste("./3_Example.outs/Figure1.2_", motif, ".boxplot.png", sep="")
## png(figfn, width=500, height=500, res=120)
## print(p)
## dev.off()

## motif <- Motifs(atac)
## data <- motif@data







