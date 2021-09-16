###
## library("rhdf5")
## library("corpcor")
library(Matrix)
## library(MASS)
## library(scales)
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

outdir <- "./1.2_DiffPeak.outs/"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F)


### Differntial analysis based reCall peaks
### last modified at 9/14/2021, by Julong wei

####################################
### 1. Generate pseudo-bulk data ###
####################################

atac <- read_rds("../1_processing/5.1_reCallPeak.outs/3_scATAC.annot.rds")

x <- granges(atac)
xrange <- ranges(x)
count <- atac@assays$ATAC@counts
anno <- data.frame(rn=rownames(count), rnz=rowSums(count),
   chr=as.character(seqnames(x)), start=start(xrange), end=end(xrange))
autosome <- as.character(1:22)
annoSel <- anno%>%dplyr::filter(rnz>0, chr%in%autosome)
Y <- count[annoSel$rn,]

meta <- atac@meta.data%>%
   mutate(treat=gsub(".*-ATAC-|_.*", "", NEW_BARCODE)) 

meta <- meta%>%
   mutate(bti=paste(MCls, treat, SNG.BEST.GUESS, sep="_"))%>%
   dplyr::select(NEW_BARCODE, bti)

dd <- meta%>%group_by(bti)%>%summarise(ncell=n(),.groups="drop")
write_rds(dd, "./1.2_DiffPeak.outs/0_ncell.rds")

##pseudo-bulk peak data
bti <- factor(meta$bti)
X <- model.matrix(~0+bti)
YtX <- Y %*% X
YtX <- as.matrix(YtX)
colnames(YtX) <- gsub("^bti", "", colnames(YtX))

###rnz>0,chr:1-22, 260,822*400
opfn <- "./1.2_DiffPeak.outs/1_YtX.comb.rds" 
write_rds(YtX, file=opfn)


### keep features with rnz>20 and autosome and conditions with ncell>20,
### 260,821 * 331
dd <- read_rds("./1.2_DiffPeak.outs/0_ncell.rds")%>%filter(ncell>20)
YtX <- read_rds("./1.2_DiffPeak.outs/1_YtX.comb.rds")
anno <- data.frame(rn=rownames(YtX), rnz=rowSums(YtX))
annoSel <- anno%>%filter(rnz>20)

YtX_sel <- YtX[annoSel$rn, dd$bti]
opfn <- "./1.2_DiffPeak.outs/1_YtX.sel.rds"
write_rds(YtX_sel, file=opfn)





#####################################################
### 2. Differential analysis for CRP, High vs Low ###
#####################################################


rm(list=ls())

#### Read data
YtX <- read_rds("./1.2_DiffPeak.outs/1_YtX.sel.rds")

bti2 <- colnames(YtX)
cvt0 <- str_split(bti2, "_", simplify=T)
cvt <- data.frame(bti=bti2, MCls=cvt0[,1], treat=cvt0[,2], sampleID=cvt0[,3])


######################
### 2.1 call DESeq ###
######################
contrast.list <- list("LPS"=c("treat", "LPS", "CTRL"),
                      "LPS-DEX"=c("treat", "LPS-DEX", "LPS"),
                      "PHA"=c("treat", "PHA", "CTRL"),
                      "PHA-DEX"=c("treat", "PHA-DEX", "PHA"))
### run DESeq
MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell", "DC")
res <- map_dfr(MCls, function(oneX){
  time0 <- Sys.time()
###
  cvt0 <- cvt%>%dplyr::filter(MCls==oneX)
  YtX0 <- YtX[,cvt0$bti]
  rnz <- rowSums(YtX0)
  YtX0 <- YtX0[rnz>0,]
##
  dds <- DESeqDataSetFromMatrix(YtX0, cvt0, ~treat)
  dds <- DESeq(dds)
  res <- contrast.list%>%map(~results(dds, contrast=.x))
  res2 <- tibble(contrast=names(res), MCls=oneX, data=map(res,tidy))%>%unnest(data)  
##
  time1 <- Sys.time()
  elapsed <- difftime(time1, time0, units="mins")
  cat(oneX, "Features:", nrow(YtX0), "Time:", elapsed, "Done\n")
  res2
})

opfn <- "./1.2_DiffPeak.outs/2.0_DESeq.results.rds"
write_rds(res, opfn)


###
### 1. MA plots
figfn <- "./1.2_DiffPeak.outs/Figure1.1_MA.png"
png(figfn, width=900, height=1200, pointsize=12, res=150)
par(mar=c(4,4,2,2),mgp=c(2,1,0))
x <- matrix(1:20, 5, 4, byrow=T)
layout(x)

fn <- "./1.2_DiffPeak.outs/2.0_DESeq.results.rds"
res <- read_rds(fn)%>%
       mutate(color=ifelse((p.adjusted<0.1)&(!is.na(p.adjusted)), T, F))

MCls <- c("Bcell",  "Monocyte", "NKcell", "Tcell", "DC")
for (oneMCl in MCls){
##1
   res2 <- res%>%dplyr::filter(MCls==oneMCl, contrast=="LPS")%>%
           dplyr::select(baseMean, estimate, color, p.value, p.adjusted)
   print(plotMA(res2[,1:3], colLine="NA", main="LPS", cex.main=1, cex.axis=0.8, cex.lab=1))
### LPS-DEX    
   res2 <- res%>%dplyr::filter(MCls==oneMCl, contrast=="LPS-DEX")%>%
          dplyr::select(baseMean, estimate, color, p.value, p.adjusted)
   print(plotMA(res2[,1:3], colLine="NA", main="LPS-DEX", cex.main=1, cex.axis=0.8, cex.lab=1))
### PHA
   res2 <- res%>%dplyr::filter(MCls==oneMCl, contrast=="PHA")%>%
           dplyr::select(baseMean, estimate, color, p.value, p.adjusted)
   print(plotMA(res2[,1:3], colLine="NA", main="PHA", cex.main=1, cex.axis=0.8, cex.lab=1))
### PHA-DEX
   res2 <- res%>%dplyr::filter(MCls==oneMCl, contrast=="PHA-DEX")%>%
           dplyr::select(baseMean, estimate, color, p.value, p.adjusted)
   print(plotMA(res2[,1:3], colLine="NA", main="PHA-DEX", cex.main=1, cex.axis=0.8, cex.lab=1))
    
   print(mtext(oneMCl, side=4, line=0.5, cex=1, col="blue")) 
}
dev.off()


###
### 2, qq plots
fn <- "./1.2_DiffPeak.outs/2.0_DESeq.results.rds"
res <- read_rds(fn)%>%as.data.frame()%>%
   drop_na(p.value)%>%
   mutate(comb=paste(MCls, contrast, sep="_"))

comb <- sort(unique(res$comb))
dfNew <- map_dfr(comb, function(ii){
  res2 <- res%>%dplyr::filter(comb==ii)
  ngene <- nrow(res2)
  res2 <- res2%>%
     arrange(p.value)%>%
     mutate(observed=-log10(p.value), expected=-log10(ppoints(ngene)))
  res2
})

###
dfNew$MCls <- gsub("DC", "z_DC", dfNew$MCls)
lab1 <- c("LPS"="LPS", "LPS-DEX"="LPS+DEX", "PHA"="PHA", "PHA-DEX"="PHA+DEX")
lab2 <- c("Bcell"="B cell", "Monocyte"= "Monocyte", "NKcell"="NK cell",
          "Tcell"="T cell", "z_DC"="DC")
p2 <- ggplot(dfNew, aes(x=expected,y=observed))+
   ggrastr::rasterise(geom_point(size=0.3, color="grey40"),dpi=300)+
   geom_abline(colour="red")+
   facet_grid(MCls~contrast, scales="free",
      labeller=labeller(contrast=lab1,MCls=lab2))+
   xlab(bquote("Expected"~ -log[10]~"("~italic(plain(P))~")"))+
   ylab(bquote("Observed"~-log[10]~"("~italic(plain(P))~")"))+
   theme_bw()+
   theme(strip.text=element_text(size=15))

figfn <- "./1.2_DiffPeak.outs/Figure1.2_qq.png"
png(figfn, width=800, height=1000, res=120)
print(p2)
dev.off()


###
### 3, canno plots


###################################
### Table of Differential peaks ###
###################################
fn <- "./1.2_DiffPeak.outs/2.0_DESeq.results.rds"
res <- read_rds(fn)%>%drop_na(p.value)

res%>%filter(p.adjusted<0.1,abs(estimate)>0.5)%>%
   group_by(MCls, contrast)%>%
   summarise(ngene=n(),.groups="drop")

x <- res%>%filter(p.adjusted<0.1)%>%group_by(MCls)%>%nest()%>%
    mutate(ngene=map(data, ~(.x)$gene))



###########################
### if enriched in DEGs ###
###########################


