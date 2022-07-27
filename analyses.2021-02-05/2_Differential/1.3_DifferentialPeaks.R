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
library(cowplot) ##,lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
library(grid)
library(gridExtra)
library(ggExtra)
library(RColorBrewer)
library(ggsci)
library(RColorBrewer)
library(viridis)
library(ggrastr)


outdir <- "./1.3_DiffPeak.outs/"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F)


### Differntial analysis based reCall peaks
### a new model considering individual effects
### last modified at 6/23/2022, by Julong wei

####################################
### 1. Generate pseudo-bulk data ###
v####################################
 
## atac <- read_rds("../1_processing/5.1_reCallPeak.outs/3_scATAC.annot.rds")

## x <- granges(atac)
## xrange <- ranges(x)
## count <- atac@assays$ATAC@counts
## anno <- data.frame(rn=rownames(count), rnz=rowSums(count),
##    chr=as.character(seqnames(x)), start=start(xrange), end=end(xrange))
## autosome <- as.character(1:22)
## annoSel <- anno%>%dplyr::filter(rnz>0, chr%in%autosome)
## Y <- count[annoSel$rn,]

## meta <- atac@meta.data%>%
##    mutate(treat=gsub(".*-ATAC-|_.*", "", NEW_BARCODE)) 

## meta <- meta%>%
##    mutate(bti=paste(MCls, treat, SNG.BEST.GUESS, sep="_"))%>%
##    dplyr::select(NEW_BARCODE, bti)

## dd <- meta%>%group_by(bti)%>%summarise(ncell=n(),.groups="drop")
## write_rds(dd, "./1.2_DiffPeak.outs/0_ncell.rds")

## ##pseudo-bulk peak data
## bti <- factor(meta$bti)
## X <- model.matrix(~0+bti)
## YtX <- Y %*% X
## YtX <- as.matrix(YtX)
## colnames(YtX) <- gsub("^bti", "", colnames(YtX))

## ###rnz>0,chr:1-22, 260,822*400
## opfn <- "./1.2_DiffPeak.outs/1_YtX.comb.rds" 
## write_rds(YtX, file=opfn)


## ### keep features with rnz>20 and autosome and conditions with ncell>20,
## ### 260,821 * 331
## dd <- read_rds("./1.2_DiffPeak.outs/0_ncell.rds")%>%filter(ncell>20)
## YtX <- read_rds("./1.2_DiffPeak.outs/1_YtX.comb.rds")
## anno <- data.frame(rn=rownames(YtX), rnz=rowSums(YtX))
## annoSel <- anno%>%filter(rnz>20)

## YtX_sel <- YtX[annoSel$rn, dd$bti]
## opfn <- "./1.2_DiffPeak.outs/1_YtX.sel.rds"
## write_rds(YtX_sel, file=opfn)





##############################################
### 2. Differential analysis for treatment ###
##############################################


## rm(list=ls())

atac <- read_rds("../1_processing/5.1_reCallPeak.outs/3_scATAC.annot.rds")

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
  ### filter outs lowly expressed peaks
  ## atac2 <- subset(atac,subset=MCls==oneX)
  ## count <- atac2@assays$ATAC@counts
  ## rpz <- rowMeans(count>0)
  ## peakSel <- rownames(count)[rpz>0]
  
  cvt0 <- cvt%>%dplyr::filter(MCls==oneX)
  YtX0 <- YtX[,cvt0$bti]
  rnz <- rowSums(YtX0)
  peakSel2 <- rownames(YtX0)[rnz>0]
  YtX0 <- YtX0[peakSel2,]
##
  dds <- DESeqDataSetFromMatrix(YtX0, cvt0, ~treat+sampleID)
  dds <- DESeq(dds)
  res <- contrast.list%>%map(~results(dds, contrast=.x))
  res2 <- tibble(contrast=names(res), MCls=oneX, data=map(res,tidy))%>%unnest(data)  
##
  time1 <- Sys.time()
  elapsed <- difftime(time1, time0, units="mins")
  cat(oneX, "Features:", nrow(YtX0), "Time:", elapsed, "Done\n")
  res2
})

opfn <- paste(outdir, "3.0_DESeq_indi.results.rds", sep="")
write_rds(res, opfn)


###
### 1. MA plots
figfn <- "./1.3_DiffPeak.outs/Figure2.1_MA.png"
png(figfn, width=900, height=1200, pointsize=12, res=150)
par(mar=c(4,4,2,2),mgp=c(2,1,0))
x <- matrix(1:20, 5, 4, byrow=T)
layout(x)

fn <- "./1.3_DiffPeak.outs/3.0_DESeq_indi.results.rds"
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
fn <- "./1.3_DiffPeak.outs/3.0_DESeq_indi.results.rds"
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

figfn <- "./1.3_DiffPeak.outs/Figure2.2_qq.png"
png(figfn, width=800, height=1000, res=120)
print(p2)
dev.off()


###
### 3, canno plots


###################################
### Table of Differential peaks ###
###################################
## fn <- "./1.3_DiffPeak.outs/3.0_DESeq_indi.results.rds"
## res <- read_rds(fn)%>%drop_na(p.value)

## res%>%filter(p.adjusted<0.1,abs(estimate)>0.5)%>%
##    group_by(MCls, contrast)%>%
##    summarise(ngene=n(),.groups="drop")

## x <- res%>%filter(p.adjusted<0.1)%>%group_by(MCls)%>%nest()%>%
##     mutate(ngene=map(data, ~(.x)$gene))

fn <- "./1.3_DiffPeak.outs/3.0_DESeq_indi.results.rds"
res2 <- read_rds(fn)%>%as.data.frame() 
                                    
res2 <- res2%>%
   mutate(is_DEG_2=ifelse(p.adjusted<0.1&abs(estimate)>0.5, 1, 0))

x <- res2%>%group_by(MCls, contrast)%>%summarise(npeaks=sum(is_DEG_2,na.rm=T), .groups="drop")%>%as.data.frame()

x%>%pivot_wider(names_from=contrast, values_from=npeaks, values_fill=0)

npeak2 <- res2%>%dplyr::filter(is_DEG_2==1)%>%dplyr::pull(gene)%>%unique()%>%length()



################
### barplots ###
################
fn <- "./1.3_DiffPeak.outs/3.0_DESeq_indi.results.rds"
res <- read_rds(fn)%>%as.data.frame()%>%mutate(direction=ifelse(estimate>0,1,0))
res2 <- res%>%dplyr::filter(p.adjusted<0.1, abs(estimate)>0.5)
sigs <- res2%>%group_by(MCls, contrast, direction)%>%
   summarise(ny=n(), .groups="drop")%>%
   mutate(ny2=ifelse(direction==0, -ny, ny))

breaks_value <- pretty(c(-12000, 15000), 8)

p <- ggplot(sigs, aes(x=MCls, y=ny2))+
   geom_bar(aes(fill=factor(MCls), alpha=factor(direction)), stat="identity")+
   scale_fill_manual(values=c("Bcell"="#4daf4a", "DC"="#828282",
      "Monocyte"="#984ea3", "NKcell"="#aa4b56", "Tcell"="#ffaa00"))+
   scale_alpha_manual(values=c("0"=0.5, "1"=1))+
   geom_hline(yintercept=0, color="grey60")+
   geom_text(aes(x=MCls, y=ny2, label=abs(ny2),
      vjust=ifelse(direction==1, -0.2, 1.2)), size=2.8)+
   scale_y_continuous("", breaks=breaks_value, limits=c(-13000,16000),
                      labels=abs(breaks_value))+
   facet_grid(~contrast,
      labeller=labeller(contrast=c("LPS"="LPS","LPS-DEX"="LPS+DEX",
                                   "PHA"="PHA", "PHA-DEX"="PHA+DEX")))+
   theme_bw()+
   theme(legend.position="none",
         axis.title.x=element_blank(),
         axis.text.x=element_text(angle=-90, hjust=0, vjust=0.5),
         plot.margin=unit(c(5.5, 15, 5.5, 5.5), "points"))

###
figfn <- "./1.3_DiffPeak.outs/Figure3.1_barplot.png"
png(filename=figfn, width=850, height=400, res=120)
print(p)
grid.text("upregulated", x=unit(0.98,"npc"), y=unit(0.7,"npc"),
                    rot=90, hjust=0.5, vjust=0.5, gp=gpar(cex=0.9))
grid.text("downregulated", x=unit(0.98,"npc"), y=unit(0.38,"npc"),
                    rot=90, hjust=0.5, vjust=0.5, gp=gpar(cex=0.9))
dev.off()
    









#####################################################################################
### compare the results between the w/o individual model and w/. individual model ###
#####################################################################################


###
### compare results between w/o individual and w/t individuals 
df1 <- read_rds("./1.2_DiffPeak.outs/2.0_DESeq.results.rds")%>%as.data.frame()
df1 <- df1%>%drop_na(estimate, p.value)%>%
    mutate(zscore=estimate/stderror, rn=paste(MCls, contrast, gene, sep="_"))%>%
    dplyr::select(MCls, contrast, rn, zscore)
###
df2 <- read_rds("./1.3_DiffPeak.outs/3.0_DESeq_indi.results.rds")%>%as.data.frame()

df2 <- df2%>%drop_na(estimate, p.value)%>%
   mutate(zscore=estimate/stderror, rn=paste(MCls, contrast, gene, sep="_"))%>%
   dplyr::select(rn, zscore)

dfcomb <- df1%>%inner_join(df2, by="rn")

p <- ggplot(dfcomb, aes(x=zscore.x, y=zscore.y))+
   rasterise(geom_point(size=0.3, colour="grey30"), dpi=300)+
   geom_abline(colour="red")+
   facet_grid(MCls~contrast,
              labeller=labeller(contrast=c("LPS"="LPS", "LPS-DEX"="LPS+DEX",
                                           "PHA"="PHA", "PHA-DEX"="PHA+DEX")))+
   xlab("z score of baseline model(w/o individaul)")+
   ylab("z score of alternative model(w/. individual)")+
   theme_bw()
 
figfn <- paste(outdir, "Figure3_scatter.plot.png", sep="")
png(figfn, width=600, height=600, res=120)
p
dev.off()


###
### barplots
###
fn <- "./1.2_DiffPeak.outs/2.0_DESeq.results.rds"
res <- read_rds(fn)%>%mutate(rn=paste(MCls, contrast, gene, sep="_"))%>%
   as.data.frame() 
                                    
res <- res%>%
   mutate(is_DEG_1=ifelse(p.adjusted<0.1&abs(estimate)>0.5, 1, 0))%>%
   dplyr::select(MCls, contrast, gene, rn, is_DEG_1)

npeak <- res%>%dplyr::filter(is_DEG_1==1)%>%dplyr::pull(gene)%>%unique()%>%length()

### 53,453 DARs

##%>%filter(qval<0.1,abs(beta)>0.5)%>%drop_na(beta)
##
 
fn <- "./1.3_DiffPeak.outs/3.0_DESeq_indi.results.rds"
res2 <- read_rds(fn)%>%mutate(rn=paste(MCls, contrast, gene, sep="_"))%>%
   as.data.frame() 
                                    
res2 <- res2%>%
   mutate(is_DEG_2=ifelse(p.adjusted<0.1&abs(estimate)>0.5, 1, 0))

x <- res2%>%group_by(MCls, contrast)%>%summarise(npeaks=sum(is_DEG_2,na.rm=T), .groups="drop")%>%as.data.frame()

x%>%pivot_wider(names_from=contrast, values_from=npeaks, values_fill=0)

npeak2 <- res2%>%dplyr::filter(is_DEG_2==1)%>%dplyr::pull(gene)%>%unique()%>%length()
###
## %>%
##    dplyr::select(rn, is_DEG_2)
res2 <- res2%>%dplyr::select(rn, is_DEG_2)

##
resComb <- res%>%full_join(res2, by="rn")
x <- resComb$is_DEG_1
y <- resComb$is_DEG_2
grp <- rep(0, length(x))
grp[x==1&y==0] <- 1
grp[x==0&y==1] <- 3
grp[x==1&y==1] <- 2
resComb$grp <- grp

##
col2 <- c("Bcell"="#4daf4a", "DC"="#828282", "Monocyte"="#984ea3", 
          "NKcell"="#aa4b56", "Tcell"="#ffaa00")
facetlab <- as_labeller(c("LPS"="LPS", "LPS-DEX"="LPS+DEX", 
                          "PHA"="PHA", "PHA-DEX"="PHA+DEX"))

###
plotdf <- resComb%>%
   dplyr::filter(grp>0)%>%
   group_by(MCls, contrast, grp)%>%summarise(ngene=n(),.groups="drop")%>%as.data.frame()

##
p <- ggplot(plotdf, aes(x=factor(MCls), y=ngene))+
    geom_bar(aes(fill=factor(MCls), alpha=factor(grp)), stat="identity")+
    scale_fill_manual("", values=col2, guide="none")+
    scale_alpha_manual("", values=c("1"=1, "2"=0.4, "3"=0.8),
       labels=c("1"="Baseline model only", "2"="Shared", "3"="Alt model only"),
       guide=guide_legend(override.aes=list(fill="#984ea3")) )+
    facet_grid(~contrast, labeller=facetlab)+   
    theme_bw()+
    theme(axis.title=element_blank(),
         ## axis.title.y=element_text(size=12),
         axis.text.x=element_text(angle=-90, hjust=0, vjust=0.5, size=10),
         axis.text.y=element_text(size=9),
         strip.text.x=element_text(size=12))   

png("./1.3_DiffPeak.outs/Figure3.3_barplot.png", width=850, height=400, res=120)
print(p)
dev.off()


