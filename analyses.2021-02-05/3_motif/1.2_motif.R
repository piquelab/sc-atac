###
library(tidyverse)
library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(SeuratObject)
library(Signac)
library(SeuratWrappers)
library(SeuratData)
library(qqman)

library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Hsapiens.1000genomes.hs37d5, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
library(ChIPseeker, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
###
library(ggplot2)
library(cowplot, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
library(grid)
library(gridExtra)
library(ggExtra)
library(RColorBrewer)
library(ggsci)
library(viridis)
library(ComplexHeatmap)
library(circlize)
theme_set(theme_grey())



outdir <- "./1.2_motif.outs/"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)

###
### enrichmenta motif analysis for DESeq2 filtering condition 0.02 pct cells

#####################################
### enrichment analysis for motif ###
#####################################


### motif
fn <- "./1_motif.outs/1_scATAC.motif.rds" 
atac <- read_rds(fn)


###
###
peaks <- read_rds("./1.3_motif.outs/pct_0.02/1_cell-type_active.peaks.rds")


### differential results
fn <- "../2_Differential/1.2_DiffPeak.outs/2.0_DESeq.results.rds"
resDP <- read_rds(fn)%>%as.data.frame()%>%mutate(condition=paste(MCls, contrast, sep="_"))

## x <- resDP%>%filter(p.adjusted<0.1, abs(estimate)>0.5)

### fisher test
cal.fisher <- function(df){
  ###  
  resfisher <- map_dfr(1:nrow(df),function(i){  
     dmat <- matrix(as.numeric(df[i,]),2,2)
     colnames(dmat) <- c("interest", "not.interest")
     rownames(dmat) <- c("in.motif", "not.motif")
     res <- fisher.test(dmat, alternative="greater")
     res2 <- data.frame(odds=res$estimate, pval.fisher=res$p.value)
     res2
  })
  resfisher
}


###
###
## count <- atac2@assays$ATAC@counts
## ncell <- rowSums(count>0)

## bg.DP <- rownames(count)[ncell>0.01*ncol(count)]

## top.DP <- unique(resDP2$gene)

## n.interest <- length(top.DP)
## n.not <- length(bg.DP)-n.interest
    
## enrich2 <- FindMotifs(
##     object=atac2,
##     features=top.DP,
##     background=bg.DP)
## ###
## df <- data.frame("interest.in.motif"=enrich2$observed,
##    "interest.not.motif"=n.interest-enrich2$observed,
##    "not.interest.in.motif"=enrich2$background-enrich2$observed)
## ##
## df <- df%>%mutate("not.interest.not.motif"=n.not-not.interest.in.motif)
## fisher <- cal.fisher(df)
## enrich2 <- cbind(enrich2,fisher)
##  enrich2$qvalue.hyper <- p.adjust(enrich2$pvalue,"BH")
##  enrich2$qvalue.fisher <- p.adjust(enrich2$pval.fisher,"BH")

## opfn <- paste(outdir, "1_enrich.motif.rds", sep="")
## write_rds(enrich2, opfn)
 
## x <- enrich2%>%filter(odds>2)
## opfn <- paste(outdir, "2_motif.txt", sep="")
## write.table(x$motif, opfn, row.names=F, quote=F, col.names=F)



    
MCls <- rep(c("Bcell", "Monocyte", "NKcell", "Tcell"), each=4)
contrast <- rep(c("LPS", "LPS-DEX", "PHA", "PHA-DEX"), times=4)
dataset <- data.frame(MCls=MCls, contrast=contrast)

###enrichment motif analysis
enriched.motif <- lapply(1:16, function(i){
###
   cell0 <- dataset[i,1]
   contrast0 <- dataset[i,2]
   cat(i, cell0, contrast0, "\n")
   ii <- grep(cell0, colnames(peaks)) 
###
   bg.DP <- peaks[peaks[,ii]==1,1] 
   res2 <- resDP%>%dplyr::filter(MCls==cell0, contrast==contrast0)
   top.DP <- res2%>%
       dplyr::filter(p.adjusted<0.1, abs(estimate)>0.5)%>%
       dplyr::pull(gene)%>%as.character()
   top.DP <- intersect(top.DP, bg.DP)
   ## 
   n.interest <- length(top.DP)
   n.not <- length(bg.DP)-n.interest
    
   if(n.interest>0){
     enrich2 <- FindMotifs(
        object=atac,
        features=top.DP,
        background=bg.DP)
     df <- data.frame("interest.in.motif"=enrich2$observed,
                     "interest.not.motif"=n.interest-enrich2$observed,
                     "not.interest.in.motif"=enrich2$background-enrich2$observed)
     df <- df%>%mutate("not.interest.not.motif"=n.not-not.interest.in.motif)
     fisher <- cal.fisher(df)
     enrich2 <- cbind(enrich2,fisher)
     enrich2$qvalue.hyper <- p.adjust(enrich2$pvalue,"BH")
     enrich2$qvalue.fisher <- p.adjust(enrich2$pval.fisher,"BH")
     enrich2$MCls <- cell0
     enrich2$contrast <- contrast0
   }else{
     enrich2 <- NA   
   }
   enrich2 
})

enriched.motif <- enriched.motif[!is.na(enriched.motif)]
enriched.motif <- do.call(rbind,enriched.motif)

opfn <- paste(outdir, "2_motif.enrich.rds", sep="")
write_rds(enriched.motif, opfn)
##
## fn <- "../2_Differential/2.2_compareRNAandATAC.outs/2_annot.ChIPseeker.rds"
## peakAnno <- read_rds(fn)%>%
##    as.data.frame()%>%
##    mutate(chr=gsub("chr", "", seqnames),
##                    peak_region=paste(chr,start,end,sep="-"))%>%
##    dplyr::select(peak_region, geneId, SYMBOL, distanceToTSS)

## ###
## resDP2 <- resDP2%>%left_join(peakAnno, by=c("gene"="peak_region"))

## top.DP <- as.character(resDP2$gene)

### test enrichment
## build a set of background peaks
## find peaks open in T cell
## atac2 <- subset(atac, subset=MCls=="Tcell")
## cell0 <- Cells(atac2)
## open.peaks <- AccessiblePeaks(atac, cells=cell0)
## ## match the overall GC content in the peak set
## meta.feature <- GetAssayData(atac, assay="ATAC", slot="meta.features")
## peak.matched <- MatchRegionStats(
##    meta.feature=meta.feature[open.peaks,],
##    query.feature=meta.feature[top.DP,],
##    n=15000)

## ### test enrichment
## enriched <- FindMotifs(
##    object=atac,
##    features=top.DP,
##    background=peak.matched)

##

###
### heatmap
## rm(list=ls())

## fn <- paste(outdir, "2_motif.enrich.rds", sep="")
## enrich <- read_rds(fn)
## enrich <- enrich%>%
##    drop_na(fold.enrichment)%>%
##    mutate(condition=paste(MCls, contrast, sep="_"))
## condition <- unique(enrich$condition)

## ### build matrix for heatmap
## motif <- enrich$motif
## names(motif) <- enrich$motif.name
## motif <- motif[!duplicated(motif)]
## ###
## mat <- lapply(condition, function(ii){
##    enrich2 <- enrich%>%dplyr::filter(condition==ii)
##    z <- enrich2$fold.enrichment
##    names(z) <- enrich2$motif
##    z[motif]
## })
## mat <- do.call(cbind, mat)
## colnames(mat) <- condition
## rownames(mat) <- motif

## ###
## b <- as.vector(mat)
## b2 <- b[b>1&b<2.66]
## breaks <- c(seq(0,1,length.out=50),
##             quantile(b2, probs=seq(0, 1, length.out=49)), 38.93) 
## col_fun <-  colorRamp2(breaks,
##    colorRampPalette(rev(brewer.pal(n=7, name="RdBu")))(100))
## column_ha <- HeatmapAnnotation(
##     celltype=c(rep("Bcell",each=3),rep(c("Monocyte", "NKcell", "Tcell"),each=4)),
##     treatment=c("LPS-DEX","PHA","PHA-DEX",
##                rep(c("LPS", "LPS-DEX", "PHA", "PHA-DEX"),times=3)),
##     col=list(celltype=c("Bcell"="#4daf4a", "Monocyte"="#984ea3",
##                     "NKcell"="#aa4b56", "Tcell"="#ffaa00"),
##              treatment=c("LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
##                         "PHA"="#a6cee3", "PHA-DEX"="#1f78b4")))
## fig <- Heatmap(mat, col=col_fun,
##    cluster_rows=T, cluster_columns=T,
##    top_annotation=column_ha,
##    heatmap_legend_param=list(title="fold.enrichment",
##       title_gp=gpar(fontsize=10),
##       labels_gp=gpar(fontsize=10)),
##    show_row_names=F, show_column_names=T,
##    column_names_gp=gpar(fontsize=10.5),
##    raster_device="png")

## figfn <- "./1_motif.outs/Figure2.1_heatmap.png"
## png(figfn, height=800, width=600, res=120)
## set.seed(0)
## fig <- draw(fig)
## dev.off()

   

## enrich%>%mutate(LFC=log2(fold.enrichment))%>%
##    dplyr::filter(qvalue.hyper<0.1, LFC>0.5)%>%
##    group_by(condition)%>%summarise(ny=n(),.groups="drop")


## ################
## ### qq plots ###
## ################

## res <- read_rds("./1_motif.outs/2_motif.enrich.rds")%>%
##    as.data.frame()%>%
##    drop_na(pvalue)%>%
##    mutate(condition=paste(MCls, contrast, sep="_"))

## condition <- sort(unique(res$condition))
## dfNew <- map_dfr(condition, function(ii){
##   res2 <- res%>%dplyr::filter(condition==ii)
##   ngene <- nrow(res2)
##   res2 <- res2%>%
##      arrange(pvalue)%>%
##       mutate(observed=-log10(pvalue), expected=-log10(ppoints(ngene)))
##   res2
## })

## lab1 <- c("LPS"="LPS", "LPS-DEX"="LPS+DEX", "PHA"="PHA", "PHA-DEX"="PHA+DEX")
## lab2 <- c("Bcell"="B cell", "Monocyte"= "Monocyte", "NKcell"="NK cell",
##    "Tcell"="T cell")
## p <- ggplot(dfNew, aes(x=expected,y=observed))+
##    ggrastr::rasterise(geom_point(size=0.3, color="grey40"),dpi=300)+
##    geom_abline(colour="red")+
##    facet_grid(MCls~contrast, scales="free",
##       labeller=labeller(contrast=lab1, MCls=lab2))+
##    xlab(bquote("Expected"~ -log[10]~"("~italic(plain(P))~")"))+
##    ylab(bquote("Observed"~-log[10]~"("~italic(plain(P))~")"))+
##    theme_bw()+
##    theme(strip.text=element_text(size=12))

## figfn <- "./1_motif.outs/Figure2.2_qq.png"
## png(figfn, width=800, height=800, res=120)
## print(p)
## dev.off()


############################################################
### motif enrichment analysis directionally, up and down ###
############################################################

rm(list=ls())

fn <- "./1_motif.outs/1_scATAC.motif.rds" 
atac <- read_rds(fn)
## motif <- Motifs(atac)
## pfm <- GetMotifData(object=motif, slot="pwm")
## motif <- SetMotifData(object=motif, slot="pwm", new.data=pfm)
## differential peaks

peaks <- read_rds("./1.3_motif.outs/pct_0.02/1_cell-type_active.peaks.rds")

fn <- "../2_Differential/1.3_DiffPeak.outs/3.0_DESeq_indi.results.rds"
resDP <- read_rds(fn)%>%as.data.frame()%>%
   mutate(direction=ifelse(estimate>0, 1, 0))

### fisher test
cal.fisher <- function(df){
  ###  
  resfisher <- map_dfr(1:nrow(df),function(i){  
     dmat <- matrix(as.numeric(df[i,]),2,2)
     colnames(dmat) <- c("interest", "not.interest")
     rownames(dmat) <- c("in.motif", "not.motif")
     res <- fisher.test(dmat, alternative="greater")
     res2 <- data.frame(odds=res$estimate, pval.fisher=res$p.value)
     res2
  })
  resfisher
}

    
MCls <- rep(c("Bcell", "Monocyte", "NKcell", "Tcell"), each=4)
contrast <- rep(c("LPS", "LPS-DEX", "PHA", "PHA-DEX"), times=4)
dataset <- data.frame(MCls=MCls, contrast=contrast)
 
###enrichment motif analysis
enriched.motif <- lapply(1:16, function(i){
###
   cell0 <- dataset[i,1]
   contrast0 <- dataset[i,2]
   cat(i, cell0, contrast0, "\n") 
   icell <- grep(cell0, colnames(peaks)) 
   bg.DP <- peaks[peaks[,icell]==1,1] 
###
   enrich <- lapply(c(0,1),function(ii){ 
      res2 <- resDP%>%
         dplyr::filter(MCls==cell0, contrast==contrast0)
      top.DP <- res2%>%
         dplyr::filter(p.adjusted<0.1, abs(estimate)>0.5, direction==ii)%>%
         dplyr::pull(gene)%>%as.character()
      top.DP <- intersect(top.DP, bg.DP)
      ##
      n.interest <- length(top.DP)
      n.not <- length(bg.DP)-n.interest
    
     if(n.interest>100){
        enrich2 <- FindMotifs(
           object=atac,
           features=top.DP,
           background=bg.DP)
        df <- data.frame("interest.in.motif"=enrich2$observed,
           "interest.not.motif"=n.interest-enrich2$observed,
           "not.interest.in.motif"=enrich2$background-enrich2$observed)
        df <- df%>%mutate("not.interest.not.motif"=n.not-not.interest.in.motif)
        fisher <- cal.fisher(df)
        enrich2 <- cbind(enrich2, fisher)
        enrich2$qvalue.hyper <- p.adjust(enrich2$pvalue,"BH")
        enrich2$qvalue.fisher <- p.adjust(enrich2$pval.fisher,"BH")
        enrich2$MCls <- cell0
        enrich2$contrast <- contrast0
        enrich2$direction <- ii
        enrich2 <- cbind(enrich2, df)
      }else{
        enrich2 <- NA   
      }
      enrich2
   })
     
   ##return results 
   if ( sum(is.na(enrich))==2){
      enrich <- NA
   }else{   
      enrich <-enrich[!is.na(enrich)]
      enrich <- do.call(rbind,enrich)
   }   
   enrich 
})

enriched.motif <- enriched.motif[!is.na(enriched.motif)]
enriched.motif <- do.call(rbind,enriched.motif)

opfn <- paste(outdir, "3_motif.enrich.direction.rds", sep="")
write_rds(enriched.motif, opfn)



###
###
### heatmap
rm(list=ls())

outdir <- "./1.2_motif.outs/"

direction2 <- c("0"="Down","1"="Up")
fn <- paste(outdir, "3_motif.enrich.direction.rds", sep="")
enrich <- read_rds(fn)
enrich <- enrich%>%
   drop_na(fold.enrichment)%>%
   mutate(dir2=direction2[as.character(direction)],
          condition=paste(MCls, contrast, dir2, sep="_"))
condition <- unique(enrich$condition)

### build matrix for heatmap
motif <- enrich$motif
names(motif) <- enrich$motif.name
motif <- motif[!duplicated(motif)]
###
mat <- lapply(condition, function(ii){
   enrich2 <- enrich%>%dplyr::filter(condition==ii)
   z <- enrich2$fold.enrichment
   names(z) <- enrich2$motif
   z[motif]
})
mat <- do.call(cbind, mat)
colnames(mat) <- condition
rownames(mat) <- motif

###
b <- as.vector(mat)
b2 <- b[b>1&b<2.69]
breaks <- c(seq(0, 1, length.out=50),
            quantile(b2, probs=seq(0, 1, length.out=49)),9.14) 
col_fun <-  colorRamp2(breaks,
   colorRampPalette(rev(brewer.pal(n=7, name="RdBu")))(100))
column_ha <- HeatmapAnnotation(
    celltype=gsub("_.*", "", condition),
    treatment=gsub(".*cell_|.*cyte_|_(D|U).*", "", condition),
    col=list(celltype=c("Bcell"="#4daf4a", "Monocyte"="#984ea3",
                    "NKcell"="#aa4b56", "Tcell"="#ffaa00"),
             treatment=c("LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
                        "PHA"="#a6cee3", "PHA-DEX"="#1f78b4")))
fig <- Heatmap(mat, col=col_fun,
   cluster_rows=T, cluster_columns=F,
   top_annotation=column_ha,
   heatmap_legend_param=list(title="fold.enrichment",
      title_gp=gpar(fontsize=10),
      labels_gp=gpar(fontsize=10)),
   show_row_names=F, show_column_names=T,
   column_names_gp=gpar(fontsize=10),
   raster_device="png")

figfn <- paste(outdir, "Figure3.1_heatmap.png", sep="")
png(figfn, height=800, width=700, res=120)
set.seed(0)
fig <- draw(fig)
dev.off()

   
 
enrich%>%mutate(LFC=log2(fold.enrichment))%>%
   dplyr::filter(qvalue.hyper<0.1, LFC>0.5)%>%
   group_by(MCls, contrast)%>%summarise(ny=n(),.groups="drop")


##############################
### barplot enriched motif ###
##############################


direction2 <- c("0"="Down","1"="Up")
fn <- paste(outdir, "3_motif.enrich.direction.rds", sep="")
enrich <- read_rds(fn)
enrich <- enrich%>%
   drop_na(fold.enrichment)%>%
   mutate(LFC=log2(fold.enrichment))

###
res <- enrich%>%dplyr::filter(qvalue.hyper<0.1, fold.enrichment>1.41)

length(unique(res$motif))

sigs <- res%>%group_by(MCls, contrast, direction)%>%
   summarise(ny=n(), .groups="drop")
###
sig4 <- sigs%>%mutate(ny2=ifelse(direction==0, -ny, ny))

breaks_value <- pretty(c(-200, 200), 5)

p <- ggplot(sig4, aes(x=MCls, y=ny2))+
   geom_bar(aes(fill=factor(MCls), alpha=factor(direction)), stat="identity")+
   scale_fill_manual(values=c("Bcell"="#4daf4a",
      "Monocyte"="#984ea3", "NKcell"="#aa4b56", "Tcell"="#ffaa00"))+
   scale_alpha_manual(values=c("0"=0.5, "1"=1))+
   geom_hline(yintercept=0, color="grey60")+
   geom_text(aes(x=MCls, y=ny2, label=abs(ny2),
      vjust=ifelse(direction==1, -0.2, 1.2)), size=3)+
   scale_y_continuous("", breaks=breaks_value, limits=c(-200,220),
                      labels=abs(breaks_value))+
   facet_grid(~contrast,
      labeller=labeller(contrast=c("LPS"="LPS","LPS-DEX"="LPS+DEX",
                                   "PHA"="PHA", "PHA-DEX"="PHA+DEX")))+
   theme_bw()+
   theme(legend.position="none",
         axis.title.x=element_blank(),
         axis.text.x=element_text(angle=-90, hjust=0, vjust=0.5))

###
figfn <- paste(outdir, "Figure3.2_barplot.png", sep="")
png(filename=figfn, width=800, height=400, pointsize=12, res=120)
print(p)
dev.off()
                      


#######################################
### dot plots show motif enrichment ###
#######################################

###
#### label
## clst <- paste(rep(c("LPS", "PHA", "LPS-DEX", "PHA-DEX"), each=4),
##    rep(c("Bcell", "Monocyte", "NKcell", "Tcell"), times=4), sep=".")
## clst <- paste(rep(clst,times=2), rep(c("Up","Down"), each=16), sep=".")
## ###
## ii <- paste(rep(c("A","B","C","D"),each=4),rep(1:4,times=4), sep="")
## clst2 <- paste(rep(c("X","Y"),each=16), rep(ii,times=2), sep=".")
## names(clst2) <- clst

## lab2 <- gsub("-", "+", clst)
## names(lab2) <- clst2 

library(ggtext, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
library(glue)

###
### prepare data for plot
fn <- paste(outdir, "3_motif.enrich.direction.rds", sep="")
enrich <- read_rds(fn)
drt2 <- c("0"="Down", "1"="Up")
enrich <- enrich%>%mutate(direction2=drt2[as.character(direction)],
   cluster=paste(contrast, MCls, direction2, sep="."))


cl <- paste( rep(c("LPS", "PHA", "LPS+DEX", "PHA+DEX"),each=4),
   rep(c("Bcell", "Monocyte", "NKcell", "Tcell"), times=4), sep=".")
cl <- paste(rep(cl,times=2), rep(c("Up","Down"), each=16), sep=".")

cluster2 <- 1:32
names(cluster2) <- cl


## ii <- paste(rep(c("A","B","C","D"),each=4),rep(1:4,times=4), sep="")
## cl2 <- paste(rep(c("X","Y"),each=16), rep(ii,times=2), sep=".")
## cluster2 <- setNames(cl2, cl)
## lab2 <- setNames(gsub("-","+",cl),cl2)
col1 <- c("LPS"="#fb9a99", "LPS+DEX"="#e31a1c", "PHA"="#a6cee3", "PHA+DEX"="#1f78b4")
col2 <- c("Bcell"="#4daf4a", "Monocyte"="#984ea3", "NKcell"="#aa4b56", "Tcell"="#ffaa00")


enrich <- enrich%>%
  mutate(cluster=gsub("-", "+", cluster),
         contrast=gsub("-", "+", contrast),
         col.contrast=col1[contrast],
         col.MCls=col2[MCls],
         ClusterValue=as.numeric(cluster2[cluster]),
         ClusterNew=glue("<i style='color:{col.contrast}'>{contrast}.<i style='color:{col.MCls}'>{MCls}.<i style='color:black'>{direction2}"),
         ClusterNew=fct_reorder(ClusterNew, ClusterValue))


##
##
## res <- read_rds("./2_motif.activities.outs/3_motif.diff.results.rds")
## topmotif2 <- res%>%dplyr::filter(qval<0.1,abs(beta)>1.41)%>%dplyr::pull(motif)%>%unique()

###
topmotif <- enrich%>%
   ## dplyr::filter(fold.enrichment>1.41, qvalue.fisher<0.1)%>%
   ## dplyr::pull(motif)%>%unique()
   group_by(cluster)%>%
   top_n(n=6, wt=fold.enrichment)%>%ungroup()%>%dplyr::pull(motif)%>%unique()

###

enrich3 <- enrich%>%
   dplyr::filter(motif%in%topmotif, fold.enrichment>1.41, qvalue.fisher<0.1)
###
p <- ggplot(enrich3, aes(x=ClusterNew, y=motif.name))+
   geom_point(aes(size=fold.enrichment, colour=qvalue.fisher))+
   scale_colour_gradient(name="FDR",
      low="blue", high="red", na.value=NA, trans="reverse", n.breaks=5,
      guide=guide_colourbar(order=1))+
   scale_size_binned("fold enrichment",
      guide=guide_bins(show.limits=T, axis=T,
          axis.show=arrow(length=unit(1.5,"mm"), ends="both"), order=2),
      n.breaks=4)+
   ## scale_x_discrete(labels=lab2)+
    
   theme_bw()+
   theme(axis.title=element_blank(),
         axis.text.x=element_markdown(angle=45, hjust=1, size=10),
         axis.text.y=element_text(size=8))
###
figfn <- paste(outdir, "Figure3.3_dotplot.png", sep="")
png(figfn, width=850, height=880, res=120)
print(p)
dev.off()

###
###
## figfn <- paste(outdir, "Figure3.3.1_dotplot.png", sep="")
## png(figfn, width=750, height=1000, res=120)
## print(p)
## dev.off()



################
### qq plots ###
################

vdrt2 <- c("0"="Down", "1"="Up")
fn <- paste(outdir, "3_motif.enrich.direction.rds", sep="")
res <- read_rds(fn)%>%
   as.data.frame()%>%
   drop_na(pvalue)%>%
   mutate(direction2=drt2[as.character(direction)],
          comb=paste(MCls, contrast, direction2, sep="_"))

comb <- sort(unique(res$comb))
dfNew <- map_dfr(comb, function(ii){
  res2 <- res%>%dplyr::filter(comb==ii)
  ngene <- nrow(res2)
  res2 <- res2%>%
     arrange(pvalue)%>%
      mutate(observed=-log10(pvalue), expected=-log10(ppoints(ngene)))
  res2
})

p <- ggplot(dfNew, aes(x=expected,y=observed))+
   ggrastr::rasterise(geom_point(size=0.3, color="grey40"),dpi=300)+
   geom_abline(colour="red")+
   facet_wrap(vars(comb), scales="free", nrow=7, ncol=4)+
   xlab(bquote("Expected"~ -log[10]~"("~italic(plain(P))~")"))+
   ylab(bquote("Observed"~-log[10]~"("~italic(plain(P))~")"))+
   theme_bw()+
   theme(axis.text=element_text(size=7),
         strip.text=element_text(size=8))

figfn <- paste(outdir, "Figure3.4_qq.png", sep="")
png(figfn, width=800, height=900, res=120)
print(p)
dev.off()




###
### example
## enrich <- read_rds("./1.2_motif.outs/3_motif.enrich.direction.rds")

ExampleGOplot <- function(cg){

### prepare data    
   ## x <- str_split(cg$GeneRatio, "/", simplify=T)
   ## GeneRatio <- as.numeric(x[,1])/as.numeric(x[,2])
   Drt2 <- c("Up"=1, "Down"=2) 
   cg <- cg%>%mutate(Direction3=Drt2[direction2],
      contrast2=paste(Direction3, contrast.x, sep="."))%>%
      mutate(contrast2=gsub("-", "+", contrast2)) 
   ## cg$size <- rep(1,nrow(cg))
   ## cg$size[GeneRatio>=0.05&GeneRatio<0.15] <- 2
   ## cg$size[GeneRatio>=0.15] <- 3 
   #
   cg <- cg%>%drop_na(odds)
   fig0 <- ggplot(cg, aes(x=contrast2, y=MCls.x))+
      geom_point(aes(size=odds, colour=qvalue.fisher))+
      scale_x_discrete(labels=c("1.LPS"="LPS.Up", "2.LPS"="LPS.Down",
         "1.LPS+DEX"="LPS+DEX.Up", "2.LPS+DEX"="LPS+DEX.Down",
         "1.PHA"="PHA.Up", "2.PHA"="PHA.Down",
         "1.PHA+DEX"="PHA+DEX.Up", "2.PHA+DEX"="PHA+DEX.Down"))+
      scale_colour_gradient(name="p.adjust",                           
         low="blue", high="red", na.value=NA, trans="reverse", n.breaks=5,
         guide=guide_colourbar(order=1))+    #"#ffa500"
      scale_size_binned("odds ratio",
         guide=guide_bins(show.limits=TRUE, axis=TRUE,
           axis.show=arrow(length=unit(1.5,"mm"), ends="both"), order=2),
         n.breaks=4)+
      theme_bw()+
      theme(axis.title=element_blank(),
         axis.text.y=element_text(size=12),
         legend.background=element_blank(),
         legend.title=element_text(size=8),
         legend.text=element_text(size=6),
         legend.key.size=grid::unit(0.6, "lines"))
   fig0
}



###
###
contrast <- c("LPS", "LPS-DEX", "PHA", "PHA-DEX")
MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
rn <- paste(rep(contrast, each=8), rep(rep(MCls, each=2), times=4),
            rep(rep(c("Down", "Up"),times=4), times=4), sep=".")
tmp <- data.frame(contrast=rep(contrast, each=8),
   MCls=rep(rep(MCls, each=2), times=4),
   direction=rep(rep(c("Down", "Up"),times=4), times=4))%>%
   mutate(rn=paste(contrast, MCls, direction, sep="."))


###
### read data
enrich <- read_rds("./1.2_motif.outs/3_motif.enrich.direction.rds")
dir2 <- c("0"="Down", "1"="Up")
enrich$direction2 <- dir2[as.character(enrich$direction)]

##
enrich2 <- enrich%>%mutate(cluster=paste(contrast, MCls, direction2, sep="."))%>%
   filter(motif.name=="NR3C1")
enrich2 <- enrich2%>%full_join(tmp, by=c("cluster"="rn"))
p1 <- ExampleGOplot(enrich2)+
    ggtitle("NR3C1")+
    theme(axis.text.x=element_text(angle=-90, size=8, hjust=0, vjust=0.5),
          plot.title=element_text(hjust=0.5))

figfn <- "./1.2_motif.outs/Figure4.1_NR3C1.png"
png(figfn, width=500, height=400, res=120)
print(p1)
dev.off()

##
enrich2 <- enrich%>%mutate(cluster=paste(contrast, MCls, direction2, sep="."))%>%
   filter(motif.name=="NR3C2")
enrich2 <- enrich2%>%full_join(tmp, by=c("cluster"="rn"))
p2 <- ExampleGOplot(enrich2)+
   ggtitle("NR3C2")+
   theme(axis.text.x=element_text(angle=-90, size=8, hjust=0, vjust=0.5),
         plot.title=element_text(hjust=0.5)) 

###
figfn <- "./1.2_motif.outs/Figure4.2_NR3C2.png"
png(figfn, width=500, height=400, res=120)
print(p2)
dev.off()
      

