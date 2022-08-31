### using cicero

## library("rhdf5")
## library("corpcor")
## library(Matrix)
## library(MASS)
## library(scales) 
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
library(cicero, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
library(monocle3)
###
library(ggplot2)
library(cowplot)
library(patchwork)
library(ggtext)
library(RColorBrewer)
library(viridis)
theme_set(theme_grey())

outdir <- "./4.2_Integrate.outs/"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)
rm(list=ls())


####
#### Final version
### cicero infer gene activity and then annotate cell-types
### Last modified, Aug 31, By Julong Wei


########################################
### infer gene activity using cicero ###
########################################

atac <- read_rds("./3_Clustering.outs/1_seurat.cluster.rds")
## create the cicero object
atac.cds <- as.cell_data_set(x=atac) ##covert to CellDataSet format
atac.cicero <- make_cicero_cds(atac.cds, reduced_coordinates=reducedDims(atac.cds)$UMAP.ATAC) ## cicero object

fn <- "/wsu/home/groups/piquelab/scHOLD/analysis.2021-03-22/1_processing/1_ATAC.outs/4_scHOLD.ATAC.clean.rds"
scHOLD <- read_rds(fn)

## Find cicero connection
genome <- seqlengths(scHOLD)
conns <- lapply(1:22, function(i){
   genome0 <- genome[i]
   cat(names(genome0), "\n")
   genome.df <- data.frame("chr"=names(genome0), "length"=genome0)
   conns <- try(run_cicero(atac.cicero, genomic_coords=genome.df),silent=T)
   if (class(conns)=="try-error") conns <- NA
   conns
})
conn2 <- conns[!is.na(conns)]
conn2 <- do.call(rbind, conn2)
write_rds(conn2, "./4.2_Integrate.outs/1.0_conns.rds")


###
## conn2 <- read_rds("./3.2_Integrate.outs/1.0_conns.rds")
conns <- conn2
ccans <- generate_ccans(conns)
## ### Add links to a seurat object
links <- ConnectionsToLinks(conns=conns, ccans=ccans)
Links(atac) <- links


### annotation
x <- Annotation(atac)
x0 <- values(x)
anno <- data.frame(chr=seqnames(x), start=start(x), end=end(x),
   strand=strand(x), feature=as.character(x0$gene_biotype),
   gene=x0$gene_id, transcript=x0$tx_id, symbol=x0$gene_name)
## anno <- anno%>%filter(strand!="*")
## make annotation
anno <- map_dfr(1:22, function(i){
   ## if (i==23) i<- "X"
   anno2 <- anno%>%filter(chr==i)
   pos <- anno2%>%filter(strand!="-")%>%
       arrange(start)%>%
       dplyr::distinct(transcript,.keep_all=T)%>%
       mutate(end=start+1)
   ###
   neg <- anno2%>%filter(strand=="-")%>%
       arrange(desc(start))%>%
       dplyr::distinct(transcript,.keep_all=T)%>%
       mutate(start=end-1)
   anno2 <- rbind(pos, neg)[,c(1:3,8)]
   names(anno2)[4] <- "gene"
   anno2
})

cds <- annotate_cds_by_site(atac.cds, anno)
###
gene.activities <- build_gene_activity_matrix(cds, conns)
atac[["ACTIVITY"]] <- CreateAssayObject(counts=gene.activities)

write_rds(atac, "./4.2_Integrate.outs/1.1_scATAC.cicero.rds")




##############################################
### integrate annotated sc-RNA and sc-ATAC ###
##############################################

atac <- read_rds("./4.2_Integrate.outs/1.1_scATAC.cicero.rds")
DefaultAssay(atac) <- "ACTIVITY"
atac <- atac%>%
   NormalizeData(assay="ACTIVITY",
   normalization.method="LogNormalize",
   scale.factor=median(atac$nCount_ACTIVITY))


###
### read sc-RNA data
## sc <- read_rds("./3_Clustering.outs/2_seurat.RNA.rds") 
## sc <- sc%>%
##    NormalizeData()%>%
##    FindVariableFeatures(nfeatures=2000)%>%
##    ScaleData()%>%
##    RunPCA(npcs=100)%>%
##    RunUMAP(dims=1:50)
### add annotation into data
## meta <- read_rds("./2_RNA.outs/3.1_scHOLD.RNAseq.annot.rds")@meta.data
## x <- meta%>%dplyr::select(barcode, predicted.celltype.l1.score, predicted.celltype.l1, predicted.celltype.l2.score, predicted.celltype.l2)
## meta <- sc@meta.data
## meta$barcode <- rownames(meta)
## metaNew <- meta%>%left_join(x, by="barcode")
## rownames(metaNew) <- metaNew$barcode
## sc <- AddMetaData(sc, metaNew)


###
### transfer cell type annotation from sc-RNAseq to sc-ATAC data
sc <- read_rds("./3_Clustering.outs/2_seurat.RNA.rds") 
anchors <- FindTransferAnchors(reference=sc, query=atac,
   features=VariableFeatures(object=sc),
   reference.assay="RNA", query.assay ="ACTIVITY", reduction="cca")
opfn <- "./4.2_Integrate.outs/2.0_anchors.rds"
write_rds(anchors, opfn)

### Using predicted.celltype.l1
pred <- TransferData(anchorset=anchors, refdata=sc$predicted.celltype.l1,
    weight.reduction=atac[["lsi"]], dims = 2:50)
pred <- pred%>%mutate(barcode=rownames(pred))
opfn <- "./4.2_Integrate.outs/2.1_pred.rds"
write_rds(pred, file=opfn)

### Using predicted.celltype.l2
pred <- TransferData(anchorset=anchors, refdata=sc$predicted.celltype.l2,
    weight.reduction=atac[["lsi"]], dims = 2:50)
pred <- pred%>%mutate(barcode=rownames(pred))
opfn <- "./4.2_Integrate.outs/2.2_pred.rds"
write_rds(pred, file=opfn)


###
### Add annotation into atac data
## atac <- read_rds("./4.2_Integrate.outs/1.1_scATAC.cicero.rds")
x <- atac@meta.data
###
pred1 <- read_rds("./4.2_Integrate.outs/2.1_pred.rds")%>%
    dplyr::select(predicted.id, barcode)%>%
    dplyr::rename("predicted.celltype.l1"="predicted.id")
###
pred2 <- read_rds("./4.2_Integrate.outs/2.2_pred.rds")%>%
    dplyr::select(predicted.id, barcode)%>%
    dplyr::rename("predicted.celltype.l2"="predicted.id")

xNew <- x%>%
   left_join(pred1, by=c("NEW_BARCODE"="barcode"))%>%
   left_join(pred2, by=c("NEW_BARCODE"="barcode"))
rownames(xNew) <- xNew$NEW_BARCODE

atac <- AddMetaData(atac, xNew)
opfn <- "./4.2_Integrate.outs/3_scATAC.annot.rds"
write_rds(atac, opfn)

### define cell type based on cluster
atac <- read_rds("./4.2_Integrate.outs/3_scATAC.annot.rds")
cell <- c("0"="Tcell", "1"="Tcell", "2"="Monocyte", "3"="Tcell",
   "4"="Bcell", "5"="NKcell", "6"="Tcell", "7"="Tcell", "8"="Tcell",
   "9"="Tcell", "10"="Tcell", "11"="DC", "12"="Tcell",
   "13"="Tcell", "14"="Tcell")
atac$MCls <- cell[as.character(atac$seurat_clusters)]
DefaultAssay(atac) <- "ATAC"
opfn <- "./4.2_Integrate.outs/3_scATAC.annot.rds"
write_rds(atac, opfn)


##################################
### summary annottaion results ###
##################################

fn <- "./4.2_Integrate.outs/3_scATAC.annot.rds"
atac <- read_rds(fn)

p0 <- DimPlot(atac, reduction="umap.atac", label=T, raster=F)+
   theme_bw()## +
   ## ## ## guides(col=guide_legend(override.aes=list(size=2),ncol=3))+
   ## theme(legend.title=element_blank(),
   ##       legend.key.size=grid::unit(0.8,"lines"))
figfn <- "./4.2_Integrate.outs/Figure1.0_atac.cluster.png"
png(figfn, width=500, height=500, res=120)
print(p0)
dev.off()

###
p1 <- DimPlot(atac, reduction="umap.atac", group.by="predicted.celltype.l1",
   label=T,  raster=F, repel=T)+
   theme_bw()+
   theme(plot.title=element_blank(),
         legend.title=element_blank(),
         ## legend.text=element_text(size=8),
         legend.key.size=grid::unit(1,"lines"))
figfn <- "./4.2_Integrate.outs/Figure1.1_atac.pred1.png"
png(figfn, width=500, height=450, res=120)
print(p1)
dev.off()

###
p2 <- DimPlot(atac, reduction="umap.atac", group.by="predicted.celltype.l2",
   label=T, label.size=2.5, raster=F, repel=T)+
   theme_bw()+
   theme(plot.title=element_blank(),
         legend.title=element_blank(),
         legend.text=element_text(size=8),
         legend.key.size=grid::unit(0.8,"lines"))

figfn <- "./4.2_Integrate.outs/Figure1.2_atac.pred2.png"
png(figfn, width=650, height=400, res=120)
print(p2)
dev.off()


###
### heatmap
fn <- "./4.2_Integrate.outs/3_scATAC.annot.rds"
atac <- read_rds(fn)
meta <- atac@meta.data


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
figfn <- "./4.2_Integrate.outs/Figure2.1_atac.heatmap.png"
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
figfn <- "./4.2_Integrate.outs/Figure2.2_atac.heatmap.png"
png(figfn, width=600, height=600, res=120)
print(p2)
dev.off()


####
####
ref <- LoadH5Seurat("../pbmc_multimodal.h5seurat")

atac <- read_rds("./4.2_Integrate.outs/1.1_scATAC.cicero.rds")
DefaultAssay(atac) <- "ACTIVITY"


################################################
### umap of cell type-related marker feature ###
################################################

fn <- "./4.2_Integrate.outs/3_scATAC.annot.rds"
atac <- read_rds(fn)
## DefaultAssay(atac) <- "ACTIVITY"

x0 <- c("MS4A1", "MS4A7", "GNLY", "CD8A")
  ## "IL7R", "CCR7", "S100A4", "CD8A", "TNFRSF18", "ID3")
##"NKG7", "CD3D"
MCls <- c("B cells", "Monocytes", "NK cells", "T cells")
anno <- data.frame(MCls=MCls, symbol=x0)

##
#### plot
figs_ls <- lapply(1:nrow(anno), function(i){
   oneMCl <- anno$MCls[i]
   gene <- anno$symbol[i]
   fig0 <- FeaturePlot(atac, features=gene)+
      scale_color_gradient(gene, low="lightgrey", high="blue")+
      ggtitle(oneMCl)+
      theme_bw()+
      theme(axis.text=element_text(size=8),
            axis.title=element_text(size=10),
            legend.title=element_text(size=8),
            legend.key.size=grid::unit(0.4,"lines"),
            legend.text=element_text(size=6),
            plot.title=element_text(size=10, hjust=0.5))
   fig0
})

###
figfn <- "./4.2_Integrate.outs/Figure3.1_feature.png"
png(figfn, width=600, height=500, res=120)
plot_grid(figs_ls[[1]], figs_ls[[2]], figs_ls[[3]], figs_ls[[4]],
   align="hv", nrow=2, ncol=2, byrow=T)
dev.off()


###
### umap color by cell types
p0 <- DimPlot(atac, reduction="umap.atac", group.by="MCls", label=T, raster=F)+
   theme_bw()+
   scale_color_manual(values=c("Bcell"="#4daf4a", "Monocyte"="#984ea3",
                               "NKcell"="#aa4b56", "Tcell"="#ffaa00",
                               "DC"="#828282"))+
   theme(legend.position="none",
         plot.title=element_blank())

   ## ## ## guides(col=guide_legend(override.aes=list(size=2),ncol=3))+
   ## theme(legend.title=element_blank(),
   ##       legend.key.size=grid::unit(0.8,"lines"))
figfn <- "./4.2_Integrate.outs/Figure3.0_umap.atac.MCls.png"
png(figfn, width=420, height=500, res=120)
print(p0)
dev.off()


p0 <- DimPlot(atac, reduction="umap.atac", group.by="MCls", label=T, raster=F)+
   theme_bw()+
   scale_color_manual(values=c("Bcell"="#4daf4a", "Monocyte"="#984ea3",
                               "NKcell"="#aa4b56", "Tcell"="#ffaa00",
                               "DC"="#828282"))+
   xlab("UMAP_1")+ylab("UMAP_2")+ 
   theme(legend.position="none",
         plot.title=element_blank())

   ## ## ## guides(col=guide_legend(override.aes=list(size=2),ncol=3))+
   ## theme(legend.title=element_blank(),
   ##       legend.key.size=grid::unit(0.8,"lines"))
figfn <- "./4.2_Integrate.outs/Figure3.0.1_umap.atac.MCls.png"
png(figfn, width=550, height=420, res=120)
print(p0)
dev.off()



######################
### coverage plots ###
######################

fn <- "./4.2_Integrate.outs/3_scATAC.annot.rds"
atac <- read_rds(fn)
## DefaultAssay(atac) <- "ACTIVITY"

geneList <- c("MS4A1", "MS4A7", "GNLY", "CD8A", "CD14", "LYZ", "IL7R", "S100A4", "CLEC9A", "CST3", "FCER1A")

### genes
##
i <- 11
geneId <- geneList[i]    
cat(i, geneId, "\n")
p <- CoveragePlot(atac, region=geneId,
   group.by="MCls", extend.upstream=1e+03, extend.downstream=1e+03, links=F)&
   scale_fill_manual(values=c("Bcell"="#4daf4a", "Monocyte"="#984ea3",
       "NKcell"="#aa4b56", "Tcell"="#ffaa00","DC"="#828282"))&
   ## ggtitle(bquote(~italic(.(geneId))))&
   theme(legend.position="none", plot.title=element_text(hjust=0.5))

###
figfn <- paste("./4.2_Integrate.outs/Figure4.", i, "_", geneId, "_coverage.png", sep="")
png(figfn, width=580, height=400, res=100)
print(p)
dev.off()
    


##
### used for AnnotationPlot plot
FindRegion <- function(
   object,
   region,
   sep = c("-", "-"),
   assay = NULL,
   extend.upstream = 0,
   extend.downstream = 0
) {
  if (!is(object = region, class2 = "GRanges")) {
      # first try to convert to coordinates, if not lookup gene
      region <- tryCatch(
         expr = suppressWarnings(
           expr = StringToGRanges(regions = region, sep = sep)
           ),
           error = function(x) {
              region <- LookupGeneCoords(
                 object = object,
                 assay = assay,
                 gene = region
                 )
                 return(region)
            }
        )
        if (is.null(x = region)) {
           stop("Gene not found")
         }
      }
      region <- suppressWarnings(expr = Extend(
          x = region,
          upstream = extend.upstream,
          downstream = extend.downstream
       )
       )
       return(region)
 }


###
###
geneList <- c("MS4A1", "FCER1A", "LYZ", "GNLY", "IL7R")
figs_ls <- lapply(1:5, function(i){
   ###
   geneId <- geneList[i]

   ### ATAC, Chromatin accessibility 
   p1 <- CoveragePlot(atac, region=geneId, show.bulk=F, peaks=F, annotation=F, 
      group.by="MCls", extend.upstream=1e+03, extend.downstream=1e+03, links=F)&
   scale_fill_manual(values=c("Bcell"="#4daf4a", "Monocyte"="#984ea3",
       "NKcell"="#aa4b56", "Tcell"="#ffaa00","DC"="#828282"))&
   ylab("Normalized accessibility")&       
   ggtitle(bquote(~italic(.(geneId))))&
   theme(legend.position="none",
         plot.title=element_text(hjust=0.5),
         axis.title.x=element_blank(),
         axis.text.x=element_blank(),
         axis.ticks.x=element_blank(),
         axis.ticks.y=element_blank(),
         strip.text.y.left=element_blank())

  ## gene annotation
  region <- FindRegion(atac, region=geneId, assay="ATAC", extend.upstream=1e+03, extend.downstream=1e+03)
  p2 <- AnnotationPlot(atac, region=region)+
     theme(axis.title.x=element_blank(),
           axis.text.x=element_blank(),
           axis.ticks.x=element_blank(),
           axis.title.y=element_blank())

  ###  
  p <- wrap_plots(p1, p2, ncol=1, heights=c(8, 1.2))  
  ##        
  p
})

###
figfn <-"./4.2_Integrate.outs/Figure5_MCls.coverage.png"
png(figfn, width=700, height=400, res=100)
plot_grid(plotlist=figs_ls, ncol=5)
dev.off()


###
### annotation
## i <- 1
## geneId <- geneList[i]
## region <- FindRegion(atac, region=geneId, assay="ATAC", extend.upstream=1e+03, extend.downstream=1e+03)
## p2 <- AnnotationPlot(atac, region=region)+
##    theme(axis.title.x=element_blank(),
##          axis.text.x=element_blank(),
##          axis.ticks.x=element_blank())
## figfn <-"./4.2_Integrate.outs/Figure5.1_MCls.coverage.png"
## png(figfn, width=580, height=400, res=100)
## p2
## dev.off()
