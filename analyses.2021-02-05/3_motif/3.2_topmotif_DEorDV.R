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
library(cowplot)
library(ggrastr)
library(grid)
library(gridExtra)
library(ggExtra)
library(RColorBrewer)
library(ggsci)
library(RColorBrewer)
library(viridis)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
## theme_set(theme_grey())

rm(list=ls())

outdir <- "./3.2_Example.outs/"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)


###
### Fun, getPeaks
## getPeaks <- function(atac, resMotif, resDP){
##    ###
##    ### motif data   
##    motif <- Motifs(atac)
##    x <- motif@data
##    peaks <- rownames(x)
    
##    motifDf <- resMotif%>%distinct(gene, .keep_all=T)%>%
##        dplyr::select(ID=gene, motif.name=motif)%>%as.data.frame()
    
##    ## motif.ID <- ConvertMotifID(motif, motif.name)
##    peaks_motif <- map_dfr(1:nrow(motifDf), function(i){ 
##       ###
##       ID <- motifDf[i,1]
##       motif.name <- motifDf[i,2] 
       
##    ### peak containing motifs
##       peak.sel <- peaks[x[,ID]==1] ## peak containing motifs
##       resDP2 <- resDP%>%dplyr::filter(gene%in%peak.sel)
##       ##DEX
##       tmp <- resMotif%>%dplyr::filter(gene==ID, grepl("DEX", contrast), qval<0.5, abs(beta)>1.41)
##       peak.DEX <- NULL
##       if ( nrow(tmp)>0){ 
##          Bmotif <- max(tmp$beta)
##          peak.DEX <- resDP%>%
##             dplyr::filter(gene%in%peak.sel,
##             sign(estimate)==sign(Bmotif), grepl("DEX", contrast))%>%dplyr::pull(gene)%>%unique()
##       }
##       ### immune stimuli
##       tmp <- resMotif%>%dplyr::filter(gene==ID, !grepl("DEX", contrast), qval<0.5, abs(beta)>1.41)
##       peak.stimuli <- NULL 
##       if ( nrow(tmp)>0){ 
##          Bmotif <- max(tmp$beta)
##          peak.stimuli <- resDP%>%
##              dplyr::filter(gene%in%peak.sel, sign(estimate)==sign(Bmotif), !grepl("DEX", contrast))%>%
##              dplyr::pull(gene)%>%unique()
##       }
##    ##    
##        peak.final <- union(peak.DEX, peak.stimuli)
##        ##
##        peakDf <- data.frame(peaks=peak.final, motif.ID=ID, motif.name=as.character(motif.name))
##        peakDf
##     })
##     ###
##     peaks_motif
## }

### cell-type specific motifs
getPeaks <- function(atac, resMotif, resDP){
   ###
   ### motif data   
   motif <- Motifs(atac)
   x <- motif@data
   peaks <- rownames(x)
    
   motifDf <- resMotif%>%distinct(gene, .keep_all=T)%>%
       dplyr::select(ID=gene, motif.name=motif)%>%as.data.frame()


   ### 
   ### loop by cell-types  
   MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell") 
   peaks_motif <- map_dfr(MCls, function(oneMCl){
   ###    
   ### loop by motifs 
   peaks_motif2 <- map_dfr(1:nrow(motifDf), function(i){ 
      ###
      ID <- motifDf[i,1]
      motif.name <- motifDf[i,2] 
       
   ### peak containing motifs
      peak.sel <- peaks[x[,ID]==1] ## peak containing motifs
      ## resDP2 <- resDP%>%dplyr::filter(gene%in%peak.sel, MCls==oneMCl)
      ##DEX
      tmp <- resMotif%>%dplyr::filter(gene==ID, grepl("DEX", contrast), MCls==oneMCl)
      peak.DEX <- NULL
      if ( nrow(tmp)>0){ 
         Bmotif <- tmp$beta[which.max(abs(tmp$beta))]
         peak.DEX <- resDP%>%
            dplyr::filter(gene%in%peak.sel, MCls==oneMCl,
            sign(estimate)==sign(Bmotif), grepl("DEX", contrast))%>%dplyr::pull(gene)%>%unique()
      }
       
      ### immune stimuli
      tmp <- resMotif%>%dplyr::filter(gene==ID, !grepl("DEX", contrast), MCls==oneMCl)
      peak.stimuli <- NULL 
      if ( nrow(tmp)>0){ 
         Bmotif <- tmp$beta[which.max(abs(tmp$beta))]
         peak.stimuli <- resDP%>%
             dplyr::filter(gene%in%peak.sel, MCls==oneMCl,
             sign(estimate)==sign(Bmotif), !grepl("DEX", contrast))%>%
             dplyr::pull(gene)%>%unique()
      }
   ##    
       peak.final <- union(peak.DEX, peak.stimuli)
       ##
       peakDf <- data.frame(MCls=oneMCl, peaks=peak.final, motif.ID=ID, motif.name=as.character(motif.name))
       peakDf
    })
    ###
    peaks_motif2
    })   
    peaks_motif
}

    
####
#### fun-1, get gene id
getGene <- function(anno, peakSel){
   ##
   anno2 <- anno%>%dplyr::filter(peak%in%peakSel, abs(distanceToTSS)<100000) ##100000
   ### gene list in a set of peaks that contain specific motif and are differentially accessible
   gene <- unique(anno2$geneId)
   gene 
}    




###########################################################
### compare motif activties and transcriptional changes ###
###########################################################


###
### read atac data 
atac <- read_rds("../3_motif/2_motif.activities.outs/1_scATAC.motifActivities.rds")
x <- atac@meta.data%>%
    mutate(treats=gsub(".*-ATAC-|_.*", "", NEW_BARCODE))
atac <- AddMetaData(atac,x)

###
### peak annotation 
fn <- "../2_Differential/2.2_compareRNAandATAC.outs/2_annot.ChIPseeker.rds"
anno <- read_rds(fn)%>%as.data.frame()%>%
       mutate(peak=paste(gsub("chr","",seqnames), start, end, sep="-"))


###
### Differential peaks results
fn <- "../2_Differential/1.3_DiffPeak.outs/3.0_DESeq_indi.results.rds"
resDP <- read_rds(fn)%>%drop_na(p.value)%>%
    mutate(comb=paste(MCls, contrast, sep="_"))%>%
    dplyr::filter(p.adjusted<0.1, abs(estimate)>0.5)%>%as.data.frame()
 


###
### motif differential results
resMotif <- read_rds("./2_motif.activities.outs/3_motif.diff.results.rds")%>%
    mutate(comb=paste(MCls, contrast, sep="_"))%>%as.data.frame()

topmotif <- resMotif%>%dplyr::filter(qval<0.1, abs(beta)>1.41)%>%dplyr::pull(gene)%>%unique()

comb <- sort(unique(resMotif$comb))

resMotif2 <- resMotif%>%dplyr::filter(gene%in%topmotif)

### motif >> peaks



###
### Differential gene expression results

### scaip6
fn <- "/nfs/rprdata/julong/SCAIP/analyses/SCAIP-B1-6_2020.03.23/6_DEG.CelltypeNew_output/Filter2/1_DESeq.results.rds"
## fn <- "/nfs/rprdata/julong/SCAIP/analyses/SCAIP-B1-6_2020.03.23/6_DEG.CelltypeNew_output/Response_reviews/3.1_DESeq.ind.results.rds"
resDG <- read_rds(fn)%>%drop_na(p.value)%>%dplyr::filter(Batch2=="SCAIP6")%>%
    mutate(comb=paste(MCls, contrast, sep="_"))%>%as.data.frame()

### meta
fn <- "/nfs/rprdata/julong/SCAIP/analyses/SCAIP-B1-6_2020.03.23/6_DEG.CelltypeNew_output/Filter2/2_meta.rds"
DE_meta <- read_rds(fn)%>%
    mutate(comb=paste(MCls, contrast, sep="_"))




#################
### plot data ###
#################

peak_motif <- getPeaks(atac, resMotif2, resDP)

comb <- sort(unique(resMotif2$comb))
 
plotDF <- map_dfr(comb, function(ii){
   ###
   oneMCl <- gsub("_.*", "", ii) 
   res2 <- resMotif2%>%dplyr::filter(comb==ii)
   DE_gene <- DE_meta%>%dplyr::filter(qval<0.1, abs(beta)>0.5)%>%dplyr::pull(gene)%>%unique() 
   LFC.RNA <- map_dfr(res2$motif, function(mm){
       ##
       peakSel <- peak_motif%>%dplyr::filter(motif.name==mm)%>%pull(peaks)%>%unique()
       gene2 <- getGene(anno, peakSel)
       LFC <- resDG%>%dplyr::filter(comb==ii, gene%in%gene2, gene%in%DE_gene)%>%
           dplyr::pull(estimate)
       ###
       if ( length(LFC)>0){
          ## 
          tmp <- data.frame(LFC.RNA=median(LFC), ngene=length(LFC)) 
       }else{
          tmp <- data.frame(LFC.RNA=0, ngene=0) 
       }
       tmp
   })
   ###
   res2 <- cbind(res2, LFC.RNA) 
   res2
})
##
plotDF <- plotDF%>%dplyr::rename("beta.x"="beta", "beta.y"="LFC.RNA")



## df2 <- plotdf%>%group_by(MCls)%>%mutate(min.x=min(beta.x), max.x=max(beta.x))%>%ungroup()%>%
##    group_by(contrast)%>%mutate(min.y=min(beta.y), max.y=max(beta.y))%>%ungroup() 
## anno_df2 <- df2%>%
##     group_by(contrast, MCls)%>%
##     nest()%>%
##     mutate(corr=map(data, ~cor.test((.x)$beta.x, (.x)$beta.y, method="pearson")),
##           eq=map(corr,feq),
##           r2=map_dbl(corr,~(.x)$estimate),
##           xpos=map_dbl(data,~xFun(.x, a=0.7)),
##           ypos=map_dbl(data,~yFun(.x, a=1)))%>%
##    dplyr::select(-data,-corr) 



### scatter plots
MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell") 
for (i in 1:length(MCls)){

## plot data
oneMCl <- MCls[i]    
plotDF2 <- plotDF%>%dplyr::filter(MCls==oneMCl) ##, gene%in%topmotif)    
plotDF2 <- plotDF2%>%mutate(gr2=ifelse(motif%in%c("NR3C1", "NFKB1"), contrast, 0))      
annoDF <- plotDF2%>%dplyr::filter(motif%in%c("NR3C1", "NFKB1"))%>%
    mutate(motif2=paste("italic(", motif, ")", sep=""))

### correlation     
corr <- cor.test(plotDF2$beta.x, plotDF2$beta.y)
r <- round(as.numeric(corr$estimate), digits=3)
pval <- corr$p.value
    
symb <- "NS"
if (pval<0.05) symb <- "*"
if (pval<0.01) symb <- "**"        
if (pval<0.001) symb <- "***"    
## annot_expression <- "paste(italic(R), \"=\", .(r),  .(symb))"
##
annot_expression <-deparse( bquote(italic(R)==~.(r)~.(symb)))
    
###
xmin <- min(plotDF2$beta.x)
xmax <- max(plotDF2$beta.x)
xpos <- xmin+0.1*(xmax-xmin)
###
ymin <- min(plotDF2$beta.y)
ymax <- max(plotDF2$beta.y)
ypos <- ymax-0.05*(ymax-ymin)    
## eq <- bquote(italic(R)==~"0.814,"~"***")
    
p <- ggplot(plotDF2)+
    rasterise(geom_point(aes(x=beta.x, y=beta.y, colour=gr2), size=0.8),
                   dpi=300)+
    annotate("text", x=xpos, y=ypos, label=annot_expression, parse=T, color="blue")+
    ggrepel::geom_text_repel(data=annoDF,
        aes(x=beta.x, y=beta.y, label=motif2, colour=gr2), size=3, parse=T,
        box.padding=0.8, max.overlaps=Inf)+
    scale_colour_manual(values=c("0"="grey50",
                                 "LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
                                 "PHA"="#a6cee3", "PHA-DEX"="#1f78b4"))+ 
    ## facet_grid(contrast~MCls, scales="fixed",
    ##    labeller=labeller(
    ##       contrast=as_labeller(c("LPS"="LPS", "LPS-DEX"="LPS+DEX",
    ##                              "PHA"="PHA", "PHA-DEX"="PHA+DEX")) ))+
    ggtitle(oneMCl)+
    scale_x_continuous("LFC on motif activities", expand=expansion(mult=0.1))+
    scale_y_continuous("LFC on motif regulated genes transcription", expand=expansion(mult=0.1))+
    theme_bw()+
    theme(legend.position="none",
          axis.title=element_text(size=12),
          plot.title=element_text(hjust=0.5))

p2 <- p+geom_smooth(data=plotDF2, aes(x=beta.x, y=beta.y),
   method="lm", formula=y~x, color="blue", size=0.5, se=F)

## slides
figfn <- paste(outdir, "Figure_DEG_1_", oneMCl, "_motif.png", sep="")
png(filename=figfn, width=380, height=380, res=100)  
print(p2)
dev.off()

cat(oneMCl, i, "\n")
    
} ###

### poster                  
## figfn <- paste(outdir, "Figure0.1.1_motif.png", sep="")
## png(filename=figfn, width=300, height=450, res=100)  
## print(p2)
## dev.off()
###
### monocytes


#######################################
### scatter plots facet by contrast ###
#######################################
 
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
   min1 <- min(dx$beta.x, na.rm=T)
   max1 <- max(dx$beta.x, na.rm=T)
   R <- max1-min1
   xpos <- min1+a*R
}
##
yFun <- function(dx,a=0.8){
   min1 <- min(dx$beta.y, na.rm=T)
   max1 <- max(dx$beta.y, na.rm=T)
   R <- max1-min1
   ypos <- min1+a*R
}


### 
col_MCls <- c("Bcell"="#4daf4a", "Monocyte"="#984ea3", "NKcell"="#aa4b56", "Tcell"="#ffaa00","DC"="#828282")
MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")

 
##
## plotDF <- plotDF%>%dplyr::rename("beta.x"="beta", "beta.y"="LFC.RNA")
 
###
for (i in c(2,4)){
## plot data
oneMCl <- MCls[i]    
plotDF2 <- plotDF%>%dplyr::filter(MCls==oneMCl) ##, gene%in%topmotif)    

anno_df2 <- plotDF2%>%
    group_by(contrast)%>%
    nest()%>%
    mutate(corr=map(data, ~cor.test((.x)$beta.x, (.x)$beta.y, method="pearson")),
          eq=map(corr,feq),
          r2=map_dbl(corr,~(.x)$estimate),
          xpos=map_dbl(data,~xFun(.x, a=0.3)),
          ypos= map_dbl(data,~yFun(.x, a=0.98)))%>%
   dplyr::select(-data,-corr)
    
    
p <- ggplot(plotDF2)+
    rasterise(geom_point(aes(x=beta.x, y=beta.y), size=0.8, color="grey50"),dpi=300)+
    geom_text(data=anno_df2, aes(x=xpos, y=ypos, label=eq), colour="blue", size=3, parse=T)+
    facet_wrap(~contrast, nrow=2, ncol=2, scales="free",
        labeller=as_labeller(c("LPS"="LPS", "LPS-DEX"="LPS+DEX", "PHA"="PHA", "PHA-DEX"="PHA+DEX")))+
    scale_x_continuous("LFC on motif activities", expand=expansion(mult=0.1))+
    scale_y_continuous("LFC on motif regulated genes transcription", expand=expansion(mult=0.1))+
    theme_bw()+
    theme(legend.position="none",
          strip.text=element_text(size=12))
         ## panel.background=element_rect(color=col_MCls[oneMCl], size=1.5))
          ##axis.title=element_text(size=12),
          ##plot.title=element_text(hjust=0.5))

p2 <- p+geom_smooth(data=plotDF2, aes(x=beta.x, y=beta.y),
   method="lm", formula=y~x, color="blue", size=0.5, se=F)

## slides
figfn <- paste(outdir, "Figure_DEG_2_",  oneMCl, "_motif_facet.png", sep="")
png(filename=figfn, width=580, height=580, res=120)  
print(p2)
dev.off()

cat(oneMCl, i, "\n")
    
} ###



##################
### top motifs ###
##################
 
## topmotif <- resMotif%>%dplyr::filter(qval<0.1, abs(beta)>1.41)%>%dplyr::pull(gene)%>%unique()
 
## for (i in c(2,4)){
## ## plot data
## oneMCl <- MCls[i]    
## plotDF2 <- plotDF%>%dplyr::filter(MCls==oneMCl, gene%in%topmotif)    

## anno_df2 <- plotDF2%>%
##     group_by(contrast)%>%
##     nest()%>%
##     mutate(corr=map(data, ~cor.test((.x)$beta.x, (.x)$beta.y, method="pearson")),
##           eq=map(corr,feq),
##           r2=map_dbl(corr,~(.x)$estimate),
##           xpos=map_dbl(data,~xFun(.x, a=0.3)),
##           ypos= map_dbl(data,~yFun(.x, a=0.98)))%>%
##    dplyr::select(-data,-corr)
    
    
## p <- ggplot(plotDF2)+
##     rasterise(geom_point(aes(x=beta.x, y=beta.y), size=0.8, color="grey50"),dpi=300)+
##     geom_text(data=anno_df2, aes(x=xpos, y=ypos, label=eq), colour="blue", size=3, parse=T)+
##     facet_wrap(~contrast, nrow=2, ncol=2, scales="free",
##         labeller=as_labeller(c("LPS"="LPS", "LPS-DEX"="LPS+DEX", "PHA"="PHA", "PHA-DEX"="PHA+DEX")))+
##     scale_x_continuous("LFC on motif activities", expand=expansion(mult=0.1))+
##     scale_y_continuous("LFC on motif regulated genes transcription", expand=expansion(mult=0.1))+
##     theme_bw()+
##     theme(legend.position="none",
##           strip.text=element_text(size=12))
##          ## panel.background=element_rect(color=col_MCls[oneMCl], size=1.5))
##           ##axis.title=element_text(size=12),
##           ##plot.title=element_text(hjust=0.5))

## p2 <- p+geom_smooth(data=plotDF2, aes(x=beta.x, y=beta.y),
##    method="lm", formula=y~x, color="blue", size=0.5, se=F)

## ## slides
## figfn <- paste(outdir, "Figure4.", i, "_", oneMCl, "_motif_facet.png", sep="")
## png(filename=figfn, width=580, height=580, res=120)  
## print(p2)
## dev.off()

## cat(oneMCl, i, "\n")
    
## } ###



################
#### Heatmap ###
################

### correlation between LFC on RNA and LFC on motif activities

comb <- sort(unique(plotDF$comb))
DFcorr <- plotDF%>% ##dplyr::filter(gene%in%topmotif)%>%
    group_by(comb)%>%
    summarise(rr=cor(beta.x, beta.y),
              pval=as.numeric(cor.test(beta.x, beta.y)$p.value), .groups="drop")%>%as.data.frame()

DFcorr2 <- DFcorr%>%
    mutate(MCls=gsub("_.*", "", comb), contrast=gsub(".*_", "", comb),
           is_sig=ifelse(pval<0.05, 1, NA),
           rr2=rr*is_sig)

mat <- DFcorr2%>%pivot_wider(id_cols=MCls, names_from=contrast, values_from=rr2)
## mat2 for plot with NA value
mat2 <- as.matrix(mat[,-1])
colnames(mat2) <- c("LPS", "LPS+DEX", "PHA", "PHA+DEX")
rownames(mat2) <- c("Bcell", "Monocyte", "NKcell", "Tcell")

### mat3 for adding text in plots
mat3 <- DFcorr2%>%pivot_wider(id_cols=MCls, names_from=contrast, values_from=rr)
mat3 <- as.matrix(mat3[,-1])


###
###
###
## num <- sort(as.numeric(mat_olap))
## mybreak <- c(quantile(num[1:460], probs=seq(0, 1, length.out=5)),
##              quantile(num[461:484], probs=seq(0,1,length.out=5)))

mycol <- colorRamp2(seq(0, 1, length.out=12), colorRampPalette(c("white", "red"))(12))

p <- Heatmap(mat2, name="PCC", na_col="grey", 
    col=mycol, cluster_rows=F, cluster_columns=F,
    row_names_gp=gpar(fontsize=10), column_names_gp=gpar(fontsize=10),
    heatmap_legend_param=list(at=c(0, 0.25, 0.5, 0.75, 1),
        legend_width=grid::unit(2, "cm"), legend_height=grid::unit(5, "cm"),
        title_gp=gpar(fontsize=8, font=2), labels_gp=gpar(fontsize=8)),
    cell_fun=function(j, i, x, y, width, height, fill){
       grid.text(round(mat3[i,j],digits=3), x, y, gp=gpar(fontsize=8))
    })

figfn <- paste(outdir, "Figure_DEG_3_topmotif_heatmap.png", sep="")
png(figfn, height=480, width=520, res=120)
print(p)
dev.off()




#################
### Examples  ###
#################

## peak_motif <- getPeaks(atac, resMotif2, resDP)
## DE_gene <- DE_meta%>%dplyr::filter(qval<0.1, abs(beta)>0.5)%>%dplyr::pull(gene)%>%unique()

## motifs_ls <- c("NFKB1", "NR3C1")

## for (i in 1:2){
## ##    
## mm <- motifs_ls[i]

## peakSel <- peak_motif%>%dplyr::filter(motif.name==mm)%>%pull(peaks)%>%unique()
## gene2 <- getGene(anno, peakSel)

## ##
## fig_ls <- lapply(c("Monocyte", "Tcell"), function(oneMCl){
##   ###  
##   resDG2  <- resDG%>%dplyr::filter(MCls==oneMCl, gene%in%gene2, gene%in%DE_gene)

## p <- ggplot(resDG2, aes(x=factor(contrast), y=estimate, fill=contrast))+
##    geom_violin(width=0.8)+
##    geom_boxplot(width=0.2, color="grey", outlier.shape=NA)+
##    ylab(bquote(~Log[2]~"fold changes of gene expression"))+
##    scale_x_discrete("",
##       labels=c("LPS"="LPS", "LPS-DEX"="LPS+DEX", "PHA"="PHA", "PHA-DEX"="PHA+DEX"))+
##    scale_fill_manual(
##       values=c("LPS"="#fb9a99", "LPS-DEX"="#e31a1c", "PHA"="#a6cee3", "PHA-DEX"="#1f78b4"))+
##    ggtitle(bquote(~italic(.(mm))~" target genes"~.(oneMCl)))+ 
##    theme_bw()+
##    theme(legend.position="none",
##          plot.title=element_text(hjust=0.5),
##          axis.text=element_text(size=9))
## ###
## p
## })

## ##
## figfn <- paste(outdir, "Figure2.", i, "_", mm, "_gene_violin.png", sep="")
## png(figfn, width=800, height=480, res=120)
## plot_grid(plotlist=fig_ls)
## dev.off()

## } ###


## ###
## ### Examples
 
## outdir2 <- "./3.2_Example.outs/Examples/"
## if (!file.exists(outdir2)) dir.create(outdir2, showWarnings=F, recursive=T)


## DE_gene <- DE_meta%>%dplyr::filter(qval<0.1, abs(beta)>0.5)%>%dplyr::pull(gene)%>%unique()

## mm <- "NFKB1"
## oneMCl <- "Monocyte"

## peakSel <- peak_motif%>%dplyr::filter(motif.name==mm, MCls==oneMCl)%>%pull(peaks)%>%unique()
## gene2 <- getGene(anno, peakSel)

## resDG2  <- resDG%>%dplyr::filter(MCls==oneMCl, gene%in%gene2, gene%in%DE_gene)%>%arrange(desc(abs(estimate)))
## ###

## geneID <- bitr(unique(resDG2$gene), fromType="ENSEMBL", toType="SYMBOL", OrgDb=org.Hs.eg.db)
                
## resDG2 <- resDG2%>%left_join(geneID, by=c("gene"="ENSEMBL"))


## ###
## ###

## atac2 <- subset(atac, MCls==oneMCl)

## ### ATAC, Chromatin accessibility
## p1 <- CoveragePlot(atac, region="CSF2", show.bulk=F, peaks=F, annotation=F, 
##       group.by="treats", extend.upstream=1e+03, extend.downstream=1e+03, links=F)&
##    scale_fill_manual(values=c("Bcell"="#4daf4a", "Monocyte"="#984ea3",
##        "NKcell"="#aa4b56", "Tcell"="#ffaa00","DC"="#828282"))&
##    ylab("Normalized accessibility")&       
##    ggtitle(bquote(~italic(.(geneId))))&
##    theme(legend.position="none",
##          plot.title=element_text(hjust=0.5),
##          axis.title.x=element_blank(),
##          axis.text.x=element_blank(),
##          axis.ticks.x=element_blank(),
##          axis.ticks.y=element_blank(),
##          strip.text.y.left=element_blank())

##   ## gene annotation
##   region <- FindRegion(atac, region=geneId, assay="ATAC", extend.upstream=1e+03, extend.downstream=1e+03)
##   p2 <- AnnotationPlot(atac, region=region)+
##      theme(axis.title.x=element_blank(),
##            axis.text.x=element_blank(),
##            axis.ticks.x=element_blank(),
##            axis.title.y=element_blank())

##   ###  
##   p <- wrap_plots(p1, p2, ncol=1, heights=c(8, 1.2))  




## ####################### 
## ### motif boxplots  ###
## #######################
 
## motifs <- Motifs(atac)

## mm.id <- ConvertMotifID(motifs, name=mm)

## mat <- atac2@assays$chromvar@data


## DF <- atac2@meta.data%>%dplyr::select(NEW_BARCODE, sampleID=SNG.BEST.GUESS, treats)
## DF$y <- as.numeric(mat[mm.id,])

## p <- ggplot(DF, aes(x=treats, y=y, fill=treats))+
##    geom_violin(width=0.8)+
##    geom_boxplot(width=0.2, color="grey", outlier.shape=NA)+
##    scale_fill_manual(values=c("CTRL"="#828282", "LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
##    "PHA"="#a6cee3", "PHA-DEX"="#1f78b4"))+
##    ylab(bquote(~italic(.(mm))~"activities"))+ 
##    theme_bw()+
##    theme(legend.position="none",
##          axis.title.x=element_blank()) 


## figfn <- paste(outdir2, "Figure1_", mm, "_violin.png", sep="")
## png(figfn, width=550, height=380, res=120)
## p
## dev.off()


## ###
## ###
## mm <- "NR3C1"

## mm.id <- ConvertMotifID(motifs, name=mm)
## DF$y <- as.numeric(mat[mm.id,])

## p2 <- ggplot(DF, aes(x=treats, y=y, fill=treats))+
##    geom_violin(width=0.8)+
##    geom_boxplot(width=0.2, color="grey", outlier.shape=NA)+
##    scale_fill_manual(values=c("CTRL"="#828282", "LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
##    "PHA"="#a6cee3", "PHA-DEX"="#1f78b4"))+
##    ylab(bquote(~italic(.(mm))~"activities"))+ 
##    theme_bw()+
##    theme(legend.position="none",
##          axis.title.x=element_blank()) 

## figfn <- paste(outdir2, "Figure1.2_", mm, "_violin.png", sep="")
## png(figfn, width=550, height=380, res=120)
## p2
## dev.off()




#############################################
### Differential gene variability results ###
#############################################


fn <- "/nfs/rprdata/julong/SCAIP/analyses/SCAIP-B1-6_2020.03.23/10_RNA.Variance_output/tmp9/3_phiNew.results"
resDG <- read.table(fn, header=T)%>%drop_na(pval)%>%dplyr::filter(batch=="SCAIP6")%>%
    mutate(comb=paste(MCls, contrast, sep="_"))%>%as.data.frame()

fn <- "/nfs/rprdata/julong/SCAIP/analyses/SCAIP-B1-6_2020.03.23/10_RNA.Variance_output/tmp9/3_phiNew.meta"
DE_meta <- read.table(fn, header=T)

peak_motif <- getPeaks(atac, resMotif2, resDP)


#################
### plot data ###
#################

comb <- sort(unique(resMotif$comb))
plotDF <- map_dfr(comb, function(ii){
   ###
   oneMCl <- gsub("_.*", "", ii) 
  res2 <- resMotif2%>%dplyr::filter(comb==ii) ##, gene%in%topmotif)
  DE_gene <- DE_meta%>%dplyr::filter(qval<0.1, abs(beta)>0.5)%>%dplyr::pull(gene)%>%unique() 
   LFC.RNA <- map_dfr(res2$motif, function(mm){
       ##
       peakSel <- peak_motif%>%dplyr::filter(motif.name==mm)%>%pull(peaks)%>%unique()
       gene2 <- getGene(anno, peakSel)
       LFC <- resDG%>%dplyr::filter(comb==ii, gene%in%gene2, gene%in%DE_gene)%>%
           dplyr::pull(beta)
       ###
       if ( length(LFC)>0){
          ## 
          tmp <- data.frame(LFC.RNA=median(LFC), ngene=length(LFC)) 
       }else{
          tmp <- data.frame(LFC.RNA=0, ngene=0) 
       }
       tmp
   })
   ###
   res2 <- cbind(res2, LFC.RNA) 
   res2
})
##
plotDF <- plotDF%>%dplyr::rename("beta.x"="beta", "beta.y"="LFC.RNA")



####################################################################
### scatter plots contrast togeter for each cell-type separately ###
####################################################################
 
MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
for (i in 1:length(MCls)){

## plot data
oneMCl <- MCls[i]    
plotDF2 <- plotDF%>%dplyr::filter(MCls==oneMCl)    
plotDF2 <- plotDF2%>%mutate(gr2=ifelse(motif%in%c("NR3C1", "NFKB1"), contrast, 0))      
annoDF <- plotDF2%>%dplyr::filter(motif%in%c("NR3C1", "NFKB1"))%>%
    mutate(motif2=paste("italic(", motif, ")", sep=""))

### correlation     
corr <- cor.test(plotDF2$beta.x, plotDF2$beta.y)
r <- round(as.numeric(corr$estimate), digits=3)
pval <- corr$p.value
    
symb <- "NS"
if (pval<0.05) symb <- "*"
if (pval<0.01) symb <- "**"        
if (pval<0.001) symb <- "***"    
## annot_expression <- "paste(italic(R), \"=\", .(r),  .(symb))"
##
annot_expression <-deparse( bquote(italic(R)==~.(r)~.(symb)))
    
###
xmin <- min(plotDF2$beta.x)
xmax <- max(plotDF2$beta.x)
xpos <- xmin+0.1*(xmax-xmin)
###
ymin <- min(plotDF2$beta.y)
ymax <- max(plotDF2$beta.y)
ypos <- ymax-0.05*(ymax-ymin)    
## eq <- bquote(italic(R)==~"0.814,"~"***")
    
p <- ggplot(plotDF2)+
    rasterise(geom_point(aes(x=beta.x, y=beta.y, colour=gr2), size=0.8),
                   dpi=300)+
    annotate("text", x=xpos, y=ypos, label=annot_expression, parse=T, color="blue")+
    ggrepel::geom_text_repel(data=annoDF,
        aes(x=beta.x, y=beta.y, label=motif2, colour=gr2), size=3, parse=T,
        box.padding=0.8, max.overlaps=Inf)+
    scale_colour_manual(values=c("0"="grey50",
                                 "LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
                                 "PHA"="#a6cee3", "PHA-DEX"="#1f78b4"))+ 
    ## facet_grid(contrast~MCls, scales="fixed",
    ##    labeller=labeller(
    ##       contrast=as_labeller(c("LPS"="LPS", "LPS-DEX"="LPS+DEX",
    ##                              "PHA"="PHA", "PHA-DEX"="PHA+DEX")) ))+
    ggtitle(oneMCl)+
    scale_x_continuous("LFC on motif activities", expand=expansion(mult=0.1))+
    scale_y_continuous("LFC on motif regulated genes variability", expand=expansion(mult=0.1))+
    theme_bw()+
    theme(legend.position="none",
          axis.title=element_text(size=12),
          plot.title=element_text(hjust=0.5))

p2 <- p+geom_smooth(data=plotDF2, aes(x=beta.x, y=beta.y),
   method="lm", formula=y~x, color="blue", size=0.5, se=F)

## slides
figfn <- paste(outdir, "Figure_DVG_1_",  oneMCl, "_motif.png", sep="")
png(filename=figfn, width=380, height=380, res=100)  
print(p2)
dev.off()

cat(oneMCl, i, "\n")
    
} ###



#######################################
### scatter plots facet by contrast ###
#######################################

MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
for (i in c(2,4)){
## plot data
oneMCl <- MCls[i]    
plotDF2 <- plotDF%>%dplyr::filter(MCls==oneMCl) ##, gene%in%topmotif)    

anno_df2 <- plotDF2%>%
    group_by(contrast)%>%
    nest()%>%
    mutate(corr=map(data, ~cor.test((.x)$beta.x, (.x)$beta.y, method="pearson")),
          eq=map(corr,feq),
          r2=map_dbl(corr,~(.x)$estimate),
          xpos=map_dbl(data,~xFun(.x, a=0.3)),
          ypos= map_dbl(data,~yFun(.x, a=0.98)))%>%
   dplyr::select(-data,-corr)
    
    
p <- ggplot(plotDF2)+
    rasterise(geom_point(aes(x=beta.x, y=beta.y), size=0.8, color="grey50"),dpi=300)+
    geom_text(data=anno_df2, aes(x=xpos, y=ypos, label=eq), colour="blue", size=3, parse=T)+
    facet_wrap(~contrast, nrow=2, ncol=2, scales="free",
        labeller=as_labeller(c("LPS"="LPS", "LPS-DEX"="LPS+DEX", "PHA"="PHA", "PHA-DEX"="PHA+DEX")))+
    scale_x_continuous("LFC on motif activities", expand=expansion(mult=0.1))+
    scale_y_continuous("LFC on motif regulated genes variability", expand=expansion(mult=0.1))+
    theme_bw()+
    theme(legend.position="none",
          strip.text=element_text(size=12))
         ## panel.background=element_rect(color=col_MCls[oneMCl], size=1.5))
          ##axis.title=element_text(size=12),
          ##plot.title=element_text(hjust=0.5))

p2 <- p+geom_smooth(data=plotDF2, aes(x=beta.x, y=beta.y),
   method="lm", formula=y~x, color="blue", size=0.5, se=F)

## slides
figfn <- paste(outdir, "Figure_DVG_2_", oneMCl, "_motif_facet.png", sep="")
png(filename=figfn, width=580, height=580, res=120)  
print(p2)
dev.off()

cat(oneMCl, i, "\n")
    
} ###



################
#### Heatmap ###
################

### correlation between LFC on RNA and LFC on motif activities

comb <- sort(unique(plotDF$comb))
DFcorr <- plotDF%>% ##dplyr::filter(gene%in%topmotif)%>%
    group_by(comb)%>%
    summarise(rr=cor(beta.x, beta.y),
              pval=as.numeric(cor.test(beta.x, beta.y)$p.value), .groups="drop")%>%as.data.frame()

DFcorr2 <- DFcorr%>%
    mutate(MCls=gsub("_.*", "", comb), contrast=gsub(".*_", "", comb),
           is_sig=ifelse(pval<0.05, 1, NA),
           rr2=rr*is_sig)

mat <- DFcorr2%>%pivot_wider(id_cols=MCls, names_from=contrast, values_from=rr2)

mat2 <- as.matrix(mat[,-1])
colnames(mat2) <- c("LPS", "LPS+DEX", "PHA", "PHA+DEX")
rownames(mat2) <- c("Bcell", "Monocyte", "NKcell", "Tcell")

mat3 <- DFcorr2%>%pivot_wider(id_cols=MCls, names_from=contrast, values_from=rr)
mat3 <- as.matrix(mat3[,-1])

###
###
###
## num <- sort(as.numeric(mat_olap))
## mybreak <- c(quantile(num[1:460], probs=seq(0, 1, length.out=5)),
##              quantile(num[461:484], probs=seq(0,1,length.out=5)))

mycol <- colorRamp2(seq(-1, 1,length.out=20), colorRampPalette(rev(brewer.pal(n=7, name="RdBu")))(20))

## mycol <- colorRamp2(seq(0, 1, length.out=12), colorRampPalette(c("white", "red"))(12))

p <- Heatmap(mat2, name="PCC", na_col="grey",
    col=mycol, cluster_rows=F, cluster_columns=F,
    row_names_gp=gpar(fontsize=10), column_names_gp=gpar(fontsize=10),
    heatmap_legend_param=list(at=c(-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1),
        legend_width=grid::unit(2, "cm"),
        legend_height=grid::unit(5, "cm"),
        labels_gp=gpar(fontsize=8),
        title_gp=gpar(fontsize=8, font=2)),
    cell_fun=function(j, i, x, y, width, height, fill){
       grid.text(round(mat3[i,j],digits=3), x, y, gp=gpar(fontsize=8))
    })

figfn <- paste(outdir, "Figure_DVG_3_topmotif_heatmap.png", sep="")
png(figfn, height=480, width=520, res=120)
print(p)
dev.off()


