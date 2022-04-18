##
library(Matrix)
library(tidyverse)
library(qqman)
library(qvalue)
## 
library(DESeq2)
library(biobroom)
## library(ashr)
library(GenomicRanges)
## library(Seurat)
## library(SeuratDisk)
## library(SeuratData)
## library(Signac)
### annotation required package
library(clusterProfiler)
library(org.Hs.eg.db)
library(annotables)
## library(TxDb.Hsapiens.UCSC.hg19.knownGene)
## library(EnsDb.Hsapiens.v75)
library(ChIPseeker, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")

##
library(ggplot2)
library(cowplot, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
library(grid)
library(gridExtra)
library(ggExtra)
library(RColorBrewer)
library(ggsci)
library(RColorBrewer)
library(viridis)
library(ggrastr)
theme_set(theme_grey())


outdir <- "./6_scatter_RNAandATAC.outs/"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F)


#################
### read data ###
#################

fn <- "./1.2_DiffPeak.outs/2.0_DESeq.results.rds"
resDP <- read_rds(fn)%>%drop_na(p.value)%>%
   mutate(comb=paste(MCls, contrast, sep="_"))%>%
   dplyr::rename("peak"="gene")%>%as.data.frame() 


fn <- "./2.2_compareRNAandATAC.outs/2_annot.ChIPseeker.rds"
peakAnno <- read_rds(fn)%>%as.data.frame()%>%
   mutate(peak=paste(gsub("chr", "", seqnames), start, end, sep="-")) 
peakAnno2 <- peakAnno%>%dplyr::select(peak, DTSS=distanceToTSS, geneId, SYMBOL)



### test if DEG is DARs 
fn <- "/nfs/rprdata/julong/SCAIP/analyses/SCAIP-B1-6_2020.03.23/6_DEG.CelltypeNew_output/Filter2/2_meta.rds"
res <- read_rds(fn)%>%drop_na(pval)%>%
   mutate(comb=paste(MCls, contrast, sep="_"))%>%
   as.data.frame() 
res2 <- res%>%filter(qval<0.1, abs(beta)>0.5)
DEG <- unique(res2$gene)


### test if DVGs is DARs
fn <- "/nfs/rprdata/julong/SCAIP/analyses/SCAIP-B1-6_2020.03.23/10_RNA.Variance_output/tmp9/3_phiNew.meta"
res <- read.table(fn, header=T) %>%drop_na(pval)%>%
   mutate(comb=paste(MCls, contrast, sep="_"))
res2 <- res%>%filter(qval<0.1, abs(beta)>0.5)
DVG <- unique(res2$gene)




###
## peakAnno3 <- peakAnno2%>%filter(geneId%in%DEG)

#################################################################
### scatter plot between LFC on ATAC and LFC on transcription ###
#################################################################

###
feq <- function(x){
   r <- round(as.numeric(x$estimate),digits=3)
   p <- x$p.value
   if(p<0.001) symb <- "***"
   if(p>=0.001 & p<0.01) symb <- "**"
   if (p>=0.01 & p<0.05) symb <- "*"
   if(p>0.05) symb <- "NS"
  
   eq <- bquote(italic(R)==.(r)~","~.(symb))
   eq 
}

##
xFun <- function(dx,a=0.5){
   min1 <- dx$min.x[1]
   max1 <- dx$max.x[1]
   R <- max1-min1
   xpos <- min1+a*R
}
##
yFun <- function(dx,a=0.8){
   min1 <- dx$min.y[1]
   max1 <- dx$max.y[1]
   R <- max1-min1
   ypos <- min1+a*R
}



###
### scatter plot of LFC of DEG vs DP
fn <- "/nfs/rprdata/julong/SCAIP/analyses/SCAIP-B1-6_2020.03.23/6_DEG.CelltypeNew_output/Filter2/1_DESeq.results.rds"
resDG <- read_rds(fn)%>%drop_na(p.value)%>%filter(Batch2=="SCAIP6", gene%in%DEG)%>%
    mutate(comb=paste(MCls, contrast, sep="_"))%>%
    as.data.frame()

comb <- sort(unique(resDG$comb))
    
plotdf <- map_dfr(comb, function(ii){
    ##
    x <- resDP%>%filter(comb==ii)%>%
       left_join(peakAnno2, by="peak")%>%filter(geneId%in%DEG)
    ##
    x2 <- x%>%filter(abs(DTSS)<1e+05)%>%
       group_by(geneId)%>%summarise(beta.x=median(estimate), .groups="drop")
       
    y <- resDG%>%filter(comb==ii)%>%
       dplyr::select(gene, beta.y=estimate, contrast, MCls, comb)
    
    ##
    df <- x2%>%inner_join(y, by=c("geneId"="gene"))%>%as.data.frame()
    df
})


df2 <- plotdf%>%group_by(MCls)%>%mutate(min.x=min(beta.x), max.x=max(beta.x))%>%ungroup()%>%
   group_by(contrast)%>%mutate(min.y=min(beta.y), max.y=max(beta.y))%>%ungroup() 
anno_df2 <- df2%>%
    group_by(contrast, MCls)%>%
    nest()%>%
    mutate(corr=map(data, ~cor.test((.x)$beta.x, (.x)$beta.y, method="pearson")),
          eq=map(corr,feq),
          r2=map_dbl(corr,~(.x)$estimate),
          xpos=map_dbl(data,~xFun(.x, a=0.7)),
          ypos=map_dbl(data,~yFun(.x, a=1)))%>%
   dplyr::select(-data,-corr)

##
p <- ggplot(plotdf, aes(x=beta.x,y=beta.y))+
    rasterise(geom_point(size=0.3, color="grey50"),dpi=300)+ 
    geom_text(data=anno_df2, aes(x=xpos, y=ypos, label=eq), colour="blue", size=3, parse=T)+
    facet_grid(contrast~MCls, scales="fixed",
       labeller=labeller(
          contrast=as_labeller(c("LPS"="LPS", "LPS-DEX"="LPS+DEX",
                                 "PHA"="PHA", "PHA-DEX"="PHA+DEX")) ))+
    scale_x_continuous("LFC on chromation accessibility", expand=expansion(mult=0.1))+
    scale_y_continuous("LFC on transcription abundance", expand=expansion(mult=0.1))+
    theme_bw()+
    theme(axis.title=element_text(size=12))
p2 <- p+geom_smooth(method="lm",formula=y~x, size=0.5, se=F)
                           
figfn <- paste(outdir, "Figure1.1_DEGvsATAC.png", sep="")
png(filename=figfn, width=700, height=700, res=120)  
print(p2)
dev.off()


###
### scatter plot, LFS of DVG vs DP
fn <- "/nfs/rprdata/julong/SCAIP/analyses/SCAIP-B1-6_2020.03.23/10_RNA.Variance_output/tmp9/3_phiNew.results"
resDG <- read.table(fn, header=T)%>%drop_na(pval)%>%filter(batch=="SCAIP6", gene%in%DVG)%>%
    mutate(comb=paste(MCls, contrast, sep="_"))%>%
    as.data.frame()

comb <- sort(unique(resDG$comb))
    
plotdf <- map_dfr(comb, function(ii){
    ##
    x <- resDP%>%filter(comb==ii)%>%
       left_join(peakAnno2, by="peak")%>%filter(geneId%in%DVG)
    ##
    x2 <- x%>%filter(abs(DTSS)<1e+05)%>%
       group_by(geneId)%>%summarise(beta.x=median(estimate), .groups="drop")
       
    y <- resDG%>%filter(comb==ii)%>%
       dplyr::select(gene, beta.y=beta, contrast, MCls, comb)
    
    ##
    df <- x2%>%inner_join(y, by=c("geneId"="gene"))%>%as.data.frame()
    df
})


df2 <- plotdf%>%group_by(MCls)%>%mutate(min.x=min(beta.x), max.x=max(beta.x))%>%ungroup()%>%
   group_by(contrast)%>%mutate(min.y=min(beta.y), max.y=max(beta.y))%>%ungroup() 
anno_df2 <- df2%>%
    group_by(contrast, MCls)%>%
    nest()%>%
    mutate(corr=map(data, ~cor.test((.x)$beta.x, (.x)$beta.y, method="pearson")),
          eq=map(corr,feq),
          r2=map_dbl(corr,~(.x)$estimate),
          xpos=map_dbl(data,~xFun(.x, a=0.6)),
          ypos=map_dbl(data,~yFun(.x, a=1)))%>%
   dplyr::select(-data,-corr)

p <- ggplot(plotdf, aes(x=beta.x,y=beta.y))+
    rasterise(geom_point(size=0.3, color="grey50"),dpi=300)+ 
    geom_text(data=anno_df2, aes(x=xpos, y=ypos, label=eq), colour="blue", size=3, parse=T)+
    facet_grid(contrast~MCls, scales="fixed",
       labeller=labeller(
          contrast=as_labeller(c("LPS"="LPS", "LPS-DEX"="LPS+DEX",
                                 "PHA"="PHA", "PHA-DEX"="PHA+DEX")) ))+
    scale_x_continuous("LFC on chromation accessibility", expand=expansion(mult=0.1))+
    scale_y_continuous("LFC on transcription variability", expand=expansion(mult=0.1))+
    theme_bw()+
    theme(axis.title=element_text(size=12))
p2 <- p+geom_smooth(method="lm",formula=y~x, size=0.5, se=F)
                           
figfn <- paste(outdir, "Figure1.2_DVGvsATAC.png", sep="")
png(filename=figfn, width=700, height=700, res=120)  
print(p2)
dev.off()
