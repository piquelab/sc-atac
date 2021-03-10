##
source("../LibraryPackage.R")

outdir <- "./2_RNA.outs/"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)
rm(list=ls())


if (FALSE){

###
### Read scRNA-seq data
prefix <- "/wsu/home/groups/piquelab/scHOLD/scRNAseq/counts_2021-01-05/"
countfn <- paste(prefix, "scHOLD-RNA-1/outs/filtered_feature_bc_matrix", sep="")
counts <- Read10X(countfn)
sc <- CreateSeuratObject(counts, project="scHOLD")
opfn <- "./2_RNA.outs/1_scHOLD.RNAseq.rds"
write_rds(sc, file=opfn)

### dimensional reduction
sc <- read_rds("./2_RNA.outs/1_scHOLD.RNAseq.rds")
sc <- NormalizeData(sc)
sc <- FindVariableFeatures(sc, nfeatures=3000)
sc <- ScaleData(sc)
sc <- RunPCA(sc, npcs=100)
sc <- RunUMAP(sc, dims=1:50)
sc <- FindNeighbors(sc, dims=1:50, verbose=T)
sc <- FindClusters(sc, resolution=0.15, verbose=T)
###
opfn <- "./2_RNA.outs/2_scHOLD.RNAseq.cluster.rds"
write_rds(sc, file=opfn)
}

###
####
sc <- read_rds("./2_RNA.outs/2_scHOLD.RNAseq.cluster.rds")
fig0 <- VlnPlot(sc, features=c("nFeature_RNA", "nCount_RNA"), group.by="orig.ident", ncol=2, pt.size=0)&
        theme_bw()+
        theme(axis.title=element_blank(),
              axis.text.x=element_blank(),
              plot.title=element_text(hjust=0.5, size=10), legend.position="none")
              
figfn <- "./2_RNA.outs/Figure1_vlnplot.png"
png(figfn, width=500, height=250, res=120)
print(fig0)
dev.off()
        

### summary
### cluster results
fig0 <- DimPlot(object=sc, label=T, raster=F)+NoLegend()+theme_bw()+theme(legend.position="none")
figfn <- "./2_RNA.outs/Figure2.1_cluster.png"
png(figfn, width=400, height=500, res=120)
print(fig0)
dev.off()

#fig0 <- DimPlot(object=sc, label=F, raster=T)+NoLegend()+
#        theme_bw()+
#        theme(axis.text=element_blank(), 
#              axis.title=element_blank(),
#              axis.ticks=element_blank(), legend.position="none")
#figfn <- "./outs/Figure2.1.1_cluster.png"
#png(figfn, width=200, height=200, res=120)
#print(fig0)
#dev.off()


### features
x0 <- c("MS4A1", "CD79A", "MS4A7", "CD14", "GNLY", "NKG7", "CD3D", "CD8A")
fig0 <- FeaturePlot(object=sc, features=x0, ncol=4, raster=F)&
        theme_bw()+
           theme(legend.title=element_blank(),
                 legend.key.size=grid::unit(0.5,"lines"),
                 plot.title=element_text(size=12, hjust=0.5))
png("./outs/Figure2.2_MCls.feature.png", width=950, height=500, res=100)
print(fig0)
dev.off()

###
### Find differentially expressed features (cluster biomarkers)
sc <- read_rds("./outs/1.2_scHOLD.RNAseq.rds")
sc.markers <- FindAllMarkers(sc, only.pos=T, min.pct=0.25, logfc.threshold=0.25)
top10 <- sc.markers%>%group_by(cluster)%>%top_n(n=5, wt=avg_log2FC)
fig0 <- DoHeatmap(sc, features=top10$gene, size=3)+
        guides(color=F)+
        scale_fill_gradient2(low=rev(c('#d1e5f0','#67a9cf','#2166ac')),
                             mid="white", high=rev(c('#b2182b','#ef8a62','#fddbc7')),
                             midpoint=0, guide="colourbar", aesthetics="fill")+
        theme(legend.title=element_blank(), 
              legend.key.size=grid::unit(0.8,"lines"),
              legend.text=element_text(size=8), 
              axis.text.y=element_text(size=6))
                             
###        
figfn <- "./outs/Figure2.3.2_top5.heatmap.png"
png(figfn, width=800, height=550, res=120)
print(fig0)
dev.off()

###
#top10 <- sc.markers%>%group_by(cluster)%>%top_n(n=5, wt=avg_log2FC)
#fig0 <- DoHeatmap(sc, features=top10$gene, label=F)+
#        NoLegend()+
#        scale_fill_gradient2(low=rev(c('#d1e5f0','#67a9cf','#2166ac')),
#                             mid="white", high=rev(c('#b2182b','#ef8a62','#fddbc7')),
#                             midpoint=0, guide="colourbar", aesthetics="fill")+
#        theme(axis.text.y=element_blank())
#figfn <- "./outs/Figure2.3.2.small_top5.heatmap.png"
#png(figfn, width=400, height=300, res=120)
#print(fig0)
#dev.off()

}


###


####################################
### 3. automatic annotation cell ###
####################################

if (FALSE){

### reference data

ref <- LoadH5Seurat("../pbmc_multimodal.h5seurat")
p1 <- DimPlot(object=ref, reduction="wnn.umap", group.by="celltype.l1",
                label=T, label.size=3, raster=F, repel=T)+ 
        NoLegend()+
        theme_bw()+
        theme(legend.position="none",
              plot.title=element_text(hjust=0.5))
###
p2 <- DimPlot(object=ref, reduction="wnn.umap", group.by="celltype.l2",
                label=T, label.size=2.5, raster=F, repel=T)+ 
        NoLegend()+
        theme_bw()+
        theme(legend.position="none",
              plot.title=element_text(hjust=0.5))
###
figfn <- "./2_RNA.outs/Figure3.1_ref.umap.png"
png(figfn, width=700, height=400, res=120)
print(plot_grid(p1, p2, ncol=2))
dev.off()

###
sc <- read_rds("./outs/1.2_scHOLD.RNAseq.rds")
sc <- SCTransform(sc, verbose=FALSE)
anchors <- FindTransferAnchors(reference=ref, query=sc,
           normalization.method="SCT",
           reference.reduction="spca",
           dims=1:50)
###
sc <- MapQuery(anchorset=anchors, query=sc, reference=ref,
      refdata = list(
      celltype.l1="celltype.l1",
      celltype.l2="celltype.l2",
      predicted_ADT="ADT"),
      reference.reduction="spca", 
      reduction.model="wnn.umap")
      
opfn <- "./outs/2.1_scHOLD.RNAseq.annot.rds"
write_rds(sc, file=opfn)

### summary results
p1 <- DimPlot(sc, reduction="ref.umap", group.by="predicted.celltype.l1",
              label=T, label.size=3, raster=F, repel=T)+
              NoLegend()+
              theme_bw()+
              theme(legend.position="none",
                    plot.title=element_text(hjust=0.5))
             
p2 <- DimPlot(sc, reduction="ref.umap", group.by="predicted.celltype.l2",
              label=T, label.size=3, raster=F, repel=T)+
              NoLegend()+
              theme_bw()+
              theme(legend.position="none",
                    plot.title=element_text(hjust=0.5))
figfn <- "./outs/Figure3.2_annot.CITEseq.png"
png(figfn, width=800, height=400, res=120)                                 
print(plot_grid(p1, p2, ncol=2))
dev.off()

### summary results from RNA UMAP
fn <- "./outs/2.2_scHOLD.RNAseq.annot.rds"
sc <- read_rds(file=fn)
p1 <- DimPlot(sc, reduction="umap", group.by="predicted.celltype.l1",
              label=T, label.size=3, raster=F, repel=T)+
              NoLegend()+
              theme_bw()+
              theme(legend.position="none",
                    plot.title=element_text(hjust=0.5))
             
p2 <- DimPlot(sc, reduction="umap", group.by="predicted.celltype.l2",
              label=T, label.size=2.5, raster=F, repel=T)+
              NoLegend()+
              theme_bw()+
              theme(legend.position="none",
                    plot.title=element_text(hjust=0.5))
figfn <- "./2_RNA.outs/Figure4.2_annot.CITEseq.rnaUMAP.png"
png(figfn, width=700, height=400, res=120)                                 
print(plot_grid(p1, p2, ncol=2))
dev.off()

###
#fn <- "./outs/2.2_scHOLD.RNAseq.annot.rds"
#sc <- read_rds(file=fn)
#fig0 <- DimPlot(sc, reduction="umap", group.by="predicted.celltype.l1",
#              label=T, label.size=3, raster=F, repel=T)+
#              NoLegend()+
#              theme_bw()+
#              theme(legend.position="none", 
#              plot.title=element_blank(),
#              axis.text=element_blank(), 
#              axis.title=element_blank(),
#              axis.ticks=element_blank())
#figfn <- "./outs/Figure3.2.1_annot.CITEseq.png"
#png(figfn, width=300, height=300, res=120)
#print(fig0)
#dev.off()
}

