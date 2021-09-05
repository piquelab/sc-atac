###
library(tidyverse)
## library(parallel)
library(data.table)
## library(purrr)
library(GenomicRanges)
library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(Signac)
library(SeuratWrappers)
library(SeuratObject)
library(EnsDb.Hsapiens.v75)
## library(cicero, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
## library(monocle3)
###
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(viridis)
library(ggrastr)
theme_set(theme_grey())


rm(list=ls())

outdir <- "./3_Clustering.outs/"
if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)

##################
### clean data ###
##################
### subset SNG cells
## atac <- read_rds("./1_Merge.outs/1_seurat.merge.rds")
## demux <- read_rds("/nfs/rprdata/julong/sc-atac/demux.2021-01-23/outs.summary/1_demuxlet.rds")
## meta <- atac@meta.data
## meta <- meta%>%left_join(demux, by=c("barcode"="NEW_BARCODE"))
## meta0 <- meta%>%dplyr::filter(DROPLET.TYPE=="SNG")
## rownames(meta0) <- meta0$barcode
## ##
## atac2 <- subset(atac, cells=meta0$barcode)
## atac2@meta.data <- meta0
## opfn <- "./3_Clustering.outs/1_seurat.SNG.rds"
## write_rds(atac2, opfn)


###########################
### Clustering analysis ###
###########################

rm(list=ls())

atac <- read_rds("./2_demux.outs/2_seurat.merge.SNG.rds")
                
atac2 <- atac%>%
   RunTFIDF()%>%
   FindTopFeatures(min.cutoff="q0")%>%
   RunSVD(n=100)%>%
   RunUMAP(reduction="lsi", dims=2:50,
           reduction.name="umap.atac", reduction.key="atacUMAP")
##
atac3 <- atac2%>%
   FindNeighbors(reduction="lsi", dims=2:50)%>%
   FindClusters(resolution=0.15, verbose=FALSE, algorithm=3)       

write_rds(atac3, file="./3_Clustering.outs/1_seurat.cluster.rds")

fig1 <- DepthCor(atac2)+
   theme(plot.title=element_text(size=10),
         plot.subtitle=element_text(size=10))
figfn <- "./3_Clustering.outs/Figure1.0_depthcor.png"
png(figfn, width=600, height=400, res=120)
print(fig1)
dev.off()


###############
### summary ###
###############

atac2 <- read_rds("./3_Clustering.outs/1_seurat.cluster.rds")

fig1 <- DimPlot(object=atac2, label=TRUE, raster=F)+
    ## NoLegend()+
    theme_bw()## +
    ## theme(legend.position="none")
figfn <- "./3_Clustering.outs/Figure1.1_cluster.png"
png(figfn, width=500, height=500, res=120)
print(fig1)
dev.off()

### data for umap plot
umap <- as.data.frame(atac2[["umap.atac"]]@cell.embeddings)
x <- atac2@meta.data
df2 <- data.frame(UMAP_1=umap[,1],
   UMAP_2=umap[,2],
   seurat_clusters=x$seurat_clusters,
   treat=gsub(".*-ATAC-|_.*", "", x$NEW_BARCODE),
   Batch=gsub("-.*", "", x$NEW_BARCODE))
## alpha <- c("LPS"="a", "LPS-DEX"="b", "PHA"="c", "PHA-DEX"="d", "CTRL"="e")
## atac2$treat.alpha <- alpha[atac2$treat]


### 2,
fig2 <- ggplot(df2, aes(x=UMAP_1, y=UMAP_2, colour=factor(seurat_clusters)))+
   rasterise(geom_point(size=0.1),dpi=300)+
   facet_wrap(~factor(Batch),ncol=2)+
   ## scale_colour_manual(values=col0,
   ##     guide=guide_legend(override.aes=list(size=2)))+
   guides(col=guide_legend(override.aes=list(size=2),ncol=1))+
   theme_bw()+
   theme(legend.title=element_blank(),
         legend.background=element_blank(),
         legend.box.background=element_blank(),
         legend.key.size=grid::unit(1,"lines"),
         strip.text=element_text(size=12))
###
png("./3_Clustering.outs/Figure1.2_BATCH.png", width=700, height=450, res=120)
print(fig2)
dev.off()


### 3,
alpha <- c("LPS"="a","LPS-DEX"="b", "PHA"="c", "PHA-DEX"="d", "CTRL"="e")
df2$treat.alpha <- alpha[df2$treat]
Labtreat <- c("a"="LPS", "b"="LPS+DEX", "c"="PHA", "d"="PHA+DEX", "e"="CTRL")

fig3 <- ggplot(df2, aes(x=UMAP_1, y=UMAP_2, colour=factor(seurat_clusters)))+
   rasterise(geom_point(size=0.1), dpi=300)+
   facet_wrap(~factor(treat.alpha), ncol=3, labeller=as_labeller(Labtreat))+
   guides(col=guide_legend(override.aes=list(size=2),ncol=3))+
   theme_bw()+
   theme(strip.text=element_text(size=12),
         legend.title=element_blank(),
         legend.position=c(0.85,0.25),
         legend.background=element_blank(),
         legend.box.background=element_blank(),
         legend.key.size=grid::unit(1,"lines"))

figfn <- "./3_Clustering.outs/Figure1.3_treat.png"
png(figfn, width=700, height=650, res=120)
print(fig3)
dev.off()


## atac2$Batch <- gsub("-.*", "", colnames(atac2))
## fig2 <- DimPlot(object=atac2, raster=F)+
##    facet_grid(~Batch)+ 
##    theme_bw()+
##    theme(legend.position="none",
##          plot.title=element_blank())
## figfn <- "./3_Clustering.outs/Figure1.2_BATCH.png"
## png(figfn, width=650, height=350, res=120)
## print(fig2)
## dev.off()

## col1 <- c("e"="#828282", 
##            "a"="#fb9a99", "b"="#e31a1c",
##            "c"="#a6cee3", "d"="#1f78b4")
## atac2$treat <- gsub(".*-ATAC-|_.*", "", colnames(atac2))
## alpha <- c("LPS"="a", "LPS-DEX"="b", "PHA"="c", "PHA-DEX"="d", "CTRL"="e")
## atac2$treat.alpha <- alpha[atac2$treat]

## fig3 <- DimPlot(object=atac2, group.by="treat.alpha", raster=F)+NoLegend()+
##         facet_wrap(~treat.alpha, ncol=2, 
##                    labeller=as_labeller(c("a"="LPS", "b"="LPS-DEX", "c"="PHA", "d"="PHA-DEX", "e"="CTRL")))+
##         scale_colour_manual(values=col1)+
##         theme_bw()+
##         theme(legend.position="none",
##               plot.title=element_blank())
## figfn <- "./3_clustering.outs/Figure2.3_treat.png"
## png(figfn, width=650, height=650, res=120)
## print(fig3)
## dev.off()


                
############
#### RNA ###
############
rm(list=ls())
### annotation cell type using reference data
ref <- LoadH5Seurat("../pbmc_multimodal.h5seurat")
sc <- read_rds("/nfs/rprdata/julong/SCAIP/analyses/SCAIP-B1-6_2020.03.23/5_IdenCelltype_output/1_SCAIP.spliced.rds")
sc2 <- subset(sc, subset=BATCH=="SCAIP6")
sc2 <- SCTransform(sc2, verbose=FALSE)
anchors <- FindTransferAnchors(reference=ref, query=sc2,
   normalization.method="SCT",
   reference.reduction="spca",
   dims=1:50)
###
sc3 <- MapQuery(anchorset=anchors, query=sc2, reference=ref,
   refdata = list(celltype.l1="celltype.l1",
   celltype.l2="celltype.l2",
   predicted_ADT="ADT"),
   reference.reduction="spca",
   reduction.model="wnn.umap")
x <- sc3@meta.data

###dimensional reduction 
sc2 <- subset(sc, subset=BATCH=="SCAIP6")
sc4.0 <- sc2%>%
   NormalizeData()%>%
   FindVariableFeatures(selection.method="vst", nfeatures=2000)%>%
   ScaleData()%>%
   RunPCA(npcs=100)%>%
   RunUMAP(dims=1:50)
sc4 <- sc4.0%>%
   FindNeighbors(dims=1:50)%>%
   FindClusters(resolution=0.15)
x <- sc4@meta.data
###
metaNew <- cbind(sc3@meta.data,x[,39:40])
sc4 <- AddMetaData(sc4,metaNew)
###
opfn <- "./3_Clustering.outs/2_seurat.RNA.rds"
write_rds(sc4, file=opfn)


###############
### summary ###
###############

fn <- "./3_Clustering.outs/2_seurat.RNA.rds"
sc <- read_rds(fn)

p1 <- DimPlot(object=sc, reduction="umap", label=T, raster=F)+
   theme_bw()+
   ## guides(col=guide_legend(override.aes=list(size=2),ncol=3))+
   theme(legend.title=element_blank(),
         legend.key.size=grid::unit(1,"lines"))    
figfn <- "./3_Clustering.outs/Figure2.1_RNA.cluster.png" 
png(figfn, width=420, height=400, res=120)
print(p1)
dev.off()

###
p2 <- DimPlot(sc, reduction="umap", group.by="predicted.celltype.l1",
   label=T,  raster=F, repel=T)+
   theme_bw()+
   theme(plot.title=element_blank(),
         legend.title=element_blank(),
         ## legend.text=element_text(size=8),
         legend.key.size=grid::unit(1,"lines"))
figfn <- "./3_Clustering.outs/Figure2.2_RNA.pred1.png" 
png(figfn, width=450, height=400, res=120)
print(p2)
dev.off()

###
p3 <- DimPlot(sc, reduction="umap", group.by="predicted.celltype.l2",
   label=T, label.size=2.5, raster=F, repel=T)+
   theme_bw()+
   theme(plot.title=element_blank(),
         legend.title=element_blank(),
         legend.text=element_text(size=8),
         legend.key.size=grid::unit(0.8,"lines"))

figfn <- "./3_Clustering.outs/Figure2.3_RNA.pred2.png"
png(figfn, width=650, height=400, res=120)
print(p3)
dev.off()


###
### heatmap
fn <- "./3_Clustering.outs/2_seurat.RNA.rds"
sc <- read_rds(fn)
meta <- sc@meta.data


df1 <- meta%>%
   group_by(predicted.celltype.l1, seurat_clusters)%>%
   summarize(Freq=n(),.groups="drop")%>%
   group_by(seurat_clusters)%>%
   mutate(Perc=Freq/sum(Freq))

p1 <- ggplot(df1, aes(x=seurat_clusters, y=predicted.celltype.l1, fill=Perc))+
   geom_tile()+
   scale_fill_gradient("Fraction of cells",
      low="#ffffc8", high="#7d0025", na.value=NA)+
   xlab("Cluster")+ylab("predicted.celltype.l1")+
   theme_bw()+theme(axis.text.x=element_text(hjust=0.5, size=10))

###
figfn <- "./3_Clustering.outs/Figure3.1_heatmap.png"
png(figfn, width=800, height=600, res=120)
print(p1)
dev.off()


###
df2 <- meta%>%
   group_by(predicted.celltype.l2, seurat_clusters)%>%
   summarize(Freq=n(),.groups="drop")%>%
   group_by(seurat_clusters)%>%
   mutate(Perc=Freq/sum(Freq))

p2 <- ggplot(df2, aes(x=seurat_clusters, y=predicted.celltype.l2, fill=Perc))+
   geom_tile()+
   scale_fill_gradient("Fraction of cells",
       low="#ffffc8", high="#7d0025", na.value=NA)+
   xlab("Cluster")+ylab("predicted.celltype.l2")+
   theme_bw()+theme(axis.text.x=element_text(hjust=0.5))

###
figfn <- "./3_Clustering.outs/Figure3.2_heatmap.png"
png(figfn, width=600, height=600, res=120)
print(p2)
dev.off()

#### previous UMAP and cell type annotation
## sc <- read_rds("/nfs/rprdata/julong/SCAIP/analyses/SCAIP-B1-6_2020.03.23/5_IdenCelltype_output/4_SCAIP.MCls.Harmony.rds")
## sc2 <- subset(sc, subset=BATCH=="SCAIP6")
## p1 <- DimPlot(object=sc2, reduction="umap", label=T, raster=F)+NoLegend()+
##       theme_bw()+
##       theme(legend.position="none")

## p2 <- DimPlot(sc2, reduction="umap", group.by="MCls",
##               label=T, label.size=2.5, raster=F, repel=T)+
##       NoLegend()+
##       theme_bw()+
##       theme(legend.position="none", plot.title=element_blank())
## figfn <- "./3_Clustering.outs/Figure4.1_annot.UMAP.png"
## png(figfn, width=700, height=400, res=120)
## print(plot_grid(p1, p2, ncol=2))
## dev.off()
