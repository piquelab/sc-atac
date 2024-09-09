###
library(tidyverse)
## library(parallel)
library(data.table)
## library(purrr)
library(GenomicRanges)
library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(SeuratObject)
library(Signac)
library(SeuratWrappers)
## library(cicero, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
## library(monocle3)
## library(EnsDb.Hsapiens.v75)
###
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(viridis)
library(ggrastr)
theme_set(theme_grey())

rm(list=ls())

###
###
outdir <- "./6_pub.outs/2_supp_plots/"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)



###########################################
### supplementary plots used for paper ####
###########################################

### Mar-18-2023
### By JW

atac <- read_rds("./5.1_reCallPeak.outs/3_scATAC.annot.rds")
###

meta <- atac@meta.data

meta2 <- meta%>%mutate(treats=gsub(".*-ATAC-|_.*", "", NEW_BARCODE),
                       treat2=gsub("-EtOH", "", treats),
                       EXP=gsub("_.*","",NEW_BARCODE))


dd2 <- meta2%>%group_by(SNG.BEST.GUESS, treat2)%>%
   summarise(ncell=n(), reads=mean(nCount_ATAC), ngene=mean(nFeature_ATAC),.groups="drop")



 
## xx2 <- meta2%>%group_by(EXP)%>%
##    summarise(ncell=n(), reads=mean(nCount_ATAC), ngene=mean(nFeature_ATAC),.groups="drop")

tmp <- dd2%>%
           group_by(treat2)%>%
           summarise(nind=n(), ncell=median(ncell), reads=median(reads),ngene=median(ngene),.groups="drop")


col1 <- c("CTRL"="#828282",
   "LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
   "PHA"="#a6cee3", "PHA-DEX"="#1f78b4")
lab1 <- c("CTRL"="CTRL",
   "LPS"="LPS", "LPS-DEX"="LPS+DEX",
   "PHA"="PHA", "PHA-DEX"="PHA+DEX")


fig1 <- ggplot(dd2, aes(x=treat2, y=ncell, fill=treat2))+
   geom_violin()+xlab("")+ylab("")+
   geom_boxplot(width=0.2, color="grey", outlier.shape=NA)+ 
   ggtitle("#Cells per sample")+
   scale_fill_manual(values=col1)+
   scale_x_discrete(labels=lab1)+ylim(0,3000)+
   theme_bw()+
   theme(legend.position="none",
   plot.title=element_text(hjust=0.5, size=14),
   axis.text.x=element_text(angle=-90, hjust=0, vjust=0.5, size=12),
   axis.text.y=element_text(size=12))

###
fig2 <- ggplot(dd2,aes(x=treat2, y=reads, fill=treat2))+
   geom_violin()+xlab("")+ylab("")+ylim(3000, 12000)+
   geom_boxplot(width=0.2, color="grey", outlier.shape=NA)+     
   ggtitle("#Reads per cell")+
   scale_fill_manual(values=col1)+
   scale_x_discrete(labels=lab1)+
   theme_bw()+
   theme(legend.position="none",
         plot.title=element_text(hjust=0.5, size=14),
         axis.text.x=element_text(angle=-90, hjust=0, vjust=0.5, size=12),
         axis.text.y=element_text(size=12))

###
fig3 <- ggplot(dd2,aes(x=treat2, y=ngene, fill=treat2))+
    geom_violin()+xlab("")+ylab("")+ylim(3000,9000)+
    geom_boxplot(width=0.2, color="grey", outlier.shape=NA)+     
    ggtitle("#Features per cell")+
    scale_fill_manual(values=col1)+
    scale_x_discrete(labels=lab1)+
    theme_bw()+
    theme(legend.position="none",
          plot.title=element_text(hjust=0.5, size=14),
          axis.text.x=element_text(angle=-90, hjust=0, vjust=0.5, size=12),
          axis.text.y=element_text(size=12))


## png("./2_kb2_output/Figure2.5_violin.png", width=800, height=500, res=120) 
figfn <- paste(outdir, "FigS1_1_violin.pdf", sep="")
pdf(figfn, width=10, height=5)
print(plot_grid(fig1, fig2, fig3, labels="AUTO", label_fontface="plain", label_x=0.1,  ncol=3))
dev.off()



#################################
### supplementary tables 
################################


dd2 <- meta2%>%group_by(SNG.BEST.GUESS, treat2)%>%
   summarise(ncell=n(), reads=mean(nCount_ATAC), ngene=mean(nFeature_ATAC),.groups="drop")%>%ungroup()
dd2 <- dd2%>%mutate(comb=paste(SNG.BEST.GUESS, treat2, sep="_")) 
                    
### number of cells
x2 <- meta2%>%group_by(SNG.BEST.GUESS, treat2, MCls)%>%summarize(ncell=n(), .groups="drop")%>%
    ungroup()%>%
    mutate(comb=paste(SNG.BEST.GUESS, treat2, sep="_"))

summ <- x2%>%pivot_wider(id_cols=comb, names_from=MCls, values_from=ncell, values_fill=0)
dd2 <- dd2%>%left_join(summ, by="comb")
 
dd_summ <- dd2%>%
    dplyr::select(comb, individual=SNG.BEST.GUESS, treatment=treat2, ncell, reads, npeaks=ngene,
                  Bcell, Monocyte, NKcell, Tcell, DC)

dd_summ <- dd_summ%>%mutate(reads=round(reads, 3), npeaks=round(npeaks, 3))

###
### output
opfn <- paste(outdir, "TableS1_0_summary.tsv", sep="")
write_tsv(dd_summ, opfn) 






###############################
### UMAP split by treatment ###
###############################
 
outdir <- "./6_pub.outs/2_supp_plots/"

###
fn <- "./4.2_Integrate.outs/3_scATAC.annot.rds"
atac <- read_rds(fn)

umap <- as.data.frame(atac[["umap.atac"]]@cell.embeddings)
x <- atac@meta.data
df2 <- data.frame(UMAP_1=umap[,1],
   UMAP_2=umap[,2],
   seurat_clusters=x$seurat_clusters,
   treat=gsub(".*-ATAC-|_.*", "", x$NEW_BARCODE),
   Batch=gsub("-.*", "", x$NEW_BARCODE),
   MCls=x$MCls)

treat_val <- c("LPS"=1, "LPS-DEX"=2, "PHA"=3, "PHA-DEX"=4, "CTRL"=5)

df2 <- df2%>%
    mutate(treat_value=as.numeric(treat_val[as.character(treat)]),
           treat2=fct_reorder(treat, treat_value))               
## atac2$treat.alpha <- alpha[atac2$treat]


### 2,
p2 <- ggplot(df2, aes(x=UMAP_1, y=UMAP_2, colour=MCls))+
   rasterise(geom_point(size=0.1),dpi=300)+
   facet_wrap(~factor(Batch),ncol=2)+
   scale_colour_manual(values=c("Bcell"="#4daf4a", "Monocyte"="#984ea3", "NKcell"="#aa4b56",
                                "Tcell"="#ffaa00","DC"="#828282"),
      labels=c("Bcell"="B cell", "Monocyte"="Monocyte", "NKcell"="NK cell", "Tcell"="T cell", "DC"="DC"),
      guide=guide_legend(override.aes=list(size=2)))+
   ##guides(col=guide_legend(override.aes=list(size=2),ncol=1))+
   theme_bw()+
   theme(legend.title=element_blank(),
         legend.background=element_blank(),
         legend.box.background=element_blank(),
         legend.key.size=grid::unit(1,"lines"),
         legend.text=element_text(size=12),
         axis.text=element_text(size=12),
         axis.title=element_text(size=12),         
         strip.text=element_text(size=14))

###
figfn <- paste(outdir, "FigS1_2_umap_BATCH.pdf", sep="")
pdf(figfn, width=7.2, height=3.5)
print(p2)
dev.off()

 
### 3,
Labtreat <- c("LPS"="LPS", "LPS-DEX"="LPS+DEX", "PHA"="PHA", "PHA-DEX"="PHA+DEX", "CTRL"="CTRL")
##
p3 <- ggplot(df2, aes(x=UMAP_1, y=UMAP_2, colour=factor(MCls)))+
   rasterise(geom_point(size=0.1), dpi=300)+
   facet_wrap(~factor(treat2), ncol=3, nrow=2, dir="v", labeller=as_labeller(Labtreat))+
   scale_colour_manual(
       values=c("Bcell"="#4daf4a", "Monocyte"="#984ea3", "NKcell"="#aa4b56", "Tcell"="#ffaa00","DC"="#828282"),
      labels=c("Bcell"="B cell", "Monocyte"="Monocyte", "NKcell"="NK cell", "Tcell"="T cell", "DC"="DC"),
      guide=guide_legend(override.aes=list(size=2)))+    
   ##guides(col=guide_legend(override.aes=list(size=2),ncol=3))+
   theme_bw()+
   theme(strip.text=element_text(size=14),
         legend.title=element_blank(),
         legend.position=c(0.85,0.25),
         legend.background=element_blank(),
         legend.box.background=element_blank(),
         legend.key.size=grid::unit(1,"lines"),
         legend.text=element_text(size=12),
         axis.text=element_text(size=12),
         axis.title=element_text(size=12))

figfn <- paste(outdir, "FigS1_2_umap_treat.pdf", sep="")
pdf(figfn, width=8, height=6)
print(p3)
dev.off()


###
### 4, suggetst option

umap <- as.data.frame(atac[["umap.atac"]]@cell.embeddings)
x <- atac@meta.data

df2 <- data.frame(UMAP_1=umap[,1],
   UMAP_2=umap[,2],
   seurat_clusters=x$seurat_clusters,
   treat=gsub("-", "+", gsub(".*-ATAC-|_.*", "", x$NEW_BARCODE)),
   Batch=gsub("-.*", "", x$NEW_BARCODE),
   MCls=x$MCls)

###
treat <- sort(unique(df2$treat))
plotDF <- map_dfr(treat, function(ii){
   ##
   df2 <- df2%>%mutate(treat_col=ifelse(treat==ii, treat, "bg"))%>%arrange(desc(treat_col))
   df2$treat_facet <- ii
   df2
})    


treat_val <- c("LPS"=1, "LPS+DEX"=2, "PHA"=3, "PHA+DEX"=4, "CTRL"=5)

plotDF <- plotDF%>%
    mutate(treat_value=as.numeric(treat_val[as.character(treat_facet)]),
           treat2_facet=fct_reorder(treat_facet, treat_value))  


col1 <- c("CTRL"="#828282", "LPS"="#fb9a99", "LPS+DEX"="#e31a1c",
   "PHA"="#a6cee3", "PHA+DEX"="#1f78b4", "bg"="#c7e9c0")
### "bg"="#bdbdbd") ##, 
##
p4 <- ggplot(plotDF, aes(x=UMAP_1, y=UMAP_2))+
   rasterise(geom_point(aes(colour=factor(treat_col), alpha=treat_col), size=0.4), dpi=300)+
   facet_wrap(~factor(treat2_facet), ncol=3, nrow=2, dir="v")+
   scale_colour_manual(values=col1,
      breaks=c("bg", "CTRL", "LPS", "LPS+DEX", "PHA", "PHA+DEX"),
      labels=c("bg"="Not in condition", "CTRL"="CTRL", "LPS"="LPS", "LPS+DEX"="LPS+DEX",
             "PHA"="PHA", "PHA+DEX"="PHA+DEX"),
      guide=guide_legend(override.aes=list(size=2)))+
   scale_alpha_manual(values=c("CTRL"=1, "LPS"=1, "LPS+DEX"=1, "PHA"=1, "PHA+DEX"=1, "bg"=0.4),
                      guide="none")+
   ##guides(col=guide_legend(override.aes=list(size=2),ncol=3))+
   theme_bw()+
   theme(strip.text=element_text(size=14),
         legend.title=element_blank(),
         legend.position=c(0.85, 0.25),
         legend.background=element_blank(),
         legend.box.background=element_blank(),
         legend.key.size=grid::unit(1,"lines"),
         legend.text=element_text(size=12),
         axis.text=element_text(size=12),
         axis.title=element_text(size=12))

figfn <- paste(outdir, "FigS1_2_umap_treat_suggest.pdf", sep="")
pdf(figfn, width=8, height=6)
print(p4)
dev.off()




## col1 <- c("CTRL"="#828282", "LPS"="#fb9a99", "LPS+DEX"="#e31a1c",
##    "PHA"="#a6cee3", "PHA+DEX"="#1f78b4")
## df2 <- df2%>%mutate(treat3=gsub("-", "+", treat))
## ##
## p4 <- ggplot(df2, aes(x=UMAP_1, y=UMAP_2, colour=factor(treat3)))+
##    rasterise(geom_point(size=0.1), dpi=300)+
##    ## facet_wrap(~factor(treat2), ncol=3, nrow=2, dir="v", labeller=as_labeller(Labtreat))+    
##    scale_colour_manual(values=col1,
##       guide=guide_legend(override.aes=list(size=2)))+    
##    theme_bw()+
##    theme(strip.text=element_text(size=14),
##          legend.title=element_blank(),
##          ## legend.position=c(0.85,0.25),
##          legend.background=element_blank(),
##          legend.box.background=element_blank(),
##          legend.key.size=grid::unit(1,"lines"),
##          legend.text=element_text(size=12),
##          axis.text=element_text(size=12),
##          axis.title=element_text(size=12))

## figfn <- paste(outdir, "FigS1_2_umap_color_treat.pdf", sep="")
## pdf(figfn, width=6, height=5)
## print(p4)
## dev.off()



#################################
### LSI correlated with depth ###
#################################

atac2 <- read_rds("./3_Clustering.outs/1_seurat.cluster.rds")

fig1 <- DepthCor(atac2)+
   ggtitle("Correlation between depth and LSI")+ 
   theme(plot.title=element_text(size=14, hjust=0.5),
         plot.subtitle=element_blank(),
         axis.text=element_text(size=12),
         axis.title=element_text(size=12))
 
## figfn <- "./3_Clustering.outs/Figure1.0_depthcor.png"
figfn <- paste(outdir, "FigS1_3_depthcor.pdf", sep="")
pdf(figfn, width=5, height=4)
print(fig1)
dev.off()


############################
### cell type annotation ###
############################


fn <- "./4.2_Integrate.outs/3_scATAC.annot.rds"
atac <- read_rds(fn)

x <- atac@meta.data

###
###
## mycol2 <- c("#fb8072", "#92cd2d", "#c38dc4", "#35978f", "#fc9016", "#b15928", "#bebada", "#80b1d3", "#fccde5",
##    "#8dd3c7", "#ffff99", "#b2182b", "#828282", "#762a83")

p1 <- DimPlot(atac, reduction="umap.atac", group.by="seurat_clusters",
   repel=T, pt.size=0.2, label=T, raster=T)+
   ggtitle("Seurat clusters")+ 
   theme_bw()+
   theme(legend.position="none",
         plot.title=element_text(hjust=0.5, size=14),
         axis.title=element_text(size=12),
         axis.text=element_text(size=12))


###
p2 <- DimPlot(atac, reduction="umap.atac", group.by="predicted.celltype.l1",
   repel=T, pt.size=0.2, label=T,  raster=T)+
   ggtitle("Automatic annotation (1)")+ 
   theme_bw()+ 
   theme(plot.title=element_text(hjust=0.5, size=14),
         legend.position="none",
         axis.title=element_text(size=12),
         axis.text=element_text(size=12))
         ## legend.title=element_blank(),
         ## legend.text=element_text(size=8),
         ## legend.key.size=grid::unit(1,"lines"))


### heatmap
df1 <- x%>%
   group_by(predicted.celltype.l1, seurat_clusters)%>%
   summarize(Freq=n(),.groups="drop")%>%
   group_by(seurat_clusters)%>%
   mutate(Perc=Freq/sum(Freq))%>%ungroup()

CL_value <- c("0"=1, "1"=2, "3"=3, "6"=4, "7"=5, "8"=6, "9"=7, "10"=8, "12"=9, "13"=10, "14"=11,
              "2"=12, "4"=13, "5"=14, "11"=15)

df1 <- df1%>%mutate(CL_val=as.numeric(CL_value[as.character(seurat_clusters)]),
                    Cluster2=fct_reorder(seurat_clusters, CL_val))

p3 <- ggplot(df1, aes(x=Cluster2, y=predicted.celltype.l1, fill=Perc))+
   geom_tile()+
   scale_fill_gradient(
      low="#ffffcc", high="#e31a1c", limits=c(0,1), n.breaks=5,
      guide=guide_legend(keywidth=unit(0.5, "cm"),
                         keyheight=unit(1,"cm")))+
   xlab("Cluster")+
   ggtitle("Fraction of cells")+ 
   theme_bw()+
   theme(legend.title=element_blank(),
         plot.title=element_text(hjust=0.5, size=12),
         axis.text.x=element_text(hjust=0.5, size=12),
         axis.text.y=element_text(size=12),
         axis.title.x=element_text(size=12),
         axis.title.y=element_blank())

p1 <- as_grob(p1)
p2 <- as_grob(p2)
p3 <- as_grob(p3)
###
figfn <- paste(outdir, "FigS1_4_comb_annot.pdf", sep="")
pdf(figfn, width=11, height=4)
plot_grid(p1, p2, p3, nrow=1, ncol=3, labels="AUTO", label_fontface="plain", label_x=0.1,
          align="h", axis="tb", rel_widths=c(1, 1, 1.3))
dev.off()





#########################
### Finner annotation ###
#########################


###
p2 <- DimPlot(atac, reduction="umap.atac", group.by="predicted.celltype.l2",
   label=T, label.size=2.5, raster=T, repel=T)+
   ggtitle("Automatic annotation (2)")+ 
   theme_bw()+ 
   theme(plot.title=element_text(hjust=0.5, size=14),
         legend.position="none",
         axis.title=element_text(size=12),
         axis.text=element_text(size=12))


 
### Heatmap
df1 <- x%>%
   group_by(predicted.celltype.l2, seurat_clusters)%>%
   summarize(Freq=n(),.groups="drop")%>%
   group_by(seurat_clusters)%>%
   mutate(Perc=Freq/sum(Freq))%>%ungroup()

CL_value <- c("0"=1, "1"=2, "3"=3, "6"=4, "7"=5, "8"=6, "9"=7, "10"=8, "12"=9, "13"=10, "14"=11,
              "2"=12, "4"=13, "5"=14, "11"=15)

df1 <- df1%>%mutate(CL_val=as.numeric(CL_value[as.character(seurat_clusters)]),
                    Cluster2=fct_reorder(seurat_clusters, CL_val))

p3 <- ggplot(df1, aes(x=Cluster2, y=predicted.celltype.l2, fill=Perc))+
   geom_tile()+
   scale_fill_gradient(
      low="#ffffcc", high="#e31a1c", limits=c(0,1), n.breaks=5,
      guide=guide_legend(keywidth=unit(0.5, "cm"),
                         keyheight=unit(1,"cm")))+
   xlab("Cluster")+
   ggtitle("Fraction of cells")+ 
   theme_bw()+
   theme(legend.title=element_blank(),
         plot.title=element_text(hjust=0.5, size=12),
         axis.text.x=element_text(hjust=0.5, size=12),
         axis.text.y=element_text(size=12),
         axis.title.x=element_text(size=12),
         axis.title.y=element_blank())

p1 <- as_grob(p1)
p2 <- as_grob(p2)
p3 <- as_grob(p3)
###
figfn <- paste(outdir, "FigS1_5_comb_annot.pdf", sep="")
pdf(figfn, width=12, height=4)
plot_grid(p1, p2, p3, nrow=1, ncol=3, labels="AUTO", label_fontface="plain", label_x=0.1,
          align="h", axis="tb", rel_widths=c(1, 1, 1.4))
dev.off()









