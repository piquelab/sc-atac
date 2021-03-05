##
source("../LibraryPackage.R")

outdir <- "./outs/"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)
rm(list=ls())


#################################
### 1. processing ATAC object ###
################################# 
if (FALSE){

prefix <- "/wsu/home/groups/piquelab/scHOLD/scATAC/counts.scHold.20210301/"

### Read data
### counts
countfn <- paste(prefix, "scHOLD-ATAC/outs/filtered_peak_bc_matrix.h5", sep="") 
counts <- Read10X_h5(countfn)

### objects 
fragfn <- paste(prefix, "scHOLD-ATAC/outs/fragments.tsv.gz", sep="")
chrom_assay <- CreateChromatinAssay(
   counts=counts, 
   sep=c(":", "-"), 
   genome="hg19",
   fragments=fragfn)
   
scHOLD <- CreateSeuratObject(chrom_assay, assay="ATAC")

### Add meta data
metafn <- paste(prefix, "scHOLD-ATAC/outs/singlecell.csv", sep="")
meta2 <- fread(metafn)
###                
x <- scHOLD@meta.data
x <- x%>%mutate(barcode=rownames(x))%>%left_join(meta2)
rownames(x) <- x$barcode
scHOLD <- AddMetaData(scHOLD, x)

###add the gene information to the object
annotations <- GetGRangesFromEnsDb(ensdb=EnsDb.Hsapiens.v75)
#seqlevelsStyle(annotations) <- "1000"
genome(annotations) <- "hg19"
Annotation(scHOLD) <- annotations

scHOLD <- scHOLD%>%NucleosomeSignal()%>%TSSEnrichment(fast=F)%>%NucleosomeSignal()

scHOLD$pct_reads_in_peaks <- scHOLD$peak_region_fragments/combined$passed_filters*100
combined$blacklist_ratio <- combined$blacklist_region_fragments/combined$peak_region_fragments

opfn <- "./outs/1.1_scHOLD.ATAC.rds"
write_rds(scHOLD, file=opfn)

### clustering analysis
atac <- read_rds("./outs/1.1_scHOLD.ATAC.rds")
atac <- atac%>%
        RunTFIDF()%>%
        FindTopFeatures(min.cutoff = "q0")%>%
        RunSVD(n=100)%>%
        RunUMAP(reduction="lsi", dims=2:50, reduction.name="umap.atac", reduction.key="atacUMAP_")
atac2 <- atac%>%
         FindNeighbors(reduction="lsi", dims=2:50)%>%
         FindClusters(resolution=0.3, verbose=FALSE, algorithm=3)
opfn <- "./outs/1.1_scHOLD.ATAC.cluster.rds"
write_rds(atac2,file=opfn)
         
fig0 <- DimPlot(object=atac2, label=T, raster=F)+NoLegend()+theme_bw()+theme(legend.position="none")
figfn <- "./outs/Figure1.1_atac.cluster.png"
png(figfn, width=400, height=500, res=120)
print(fig0)
dev.off()  

#atac2 <- read_rds("./outs/1.1_scHOLD.ATAC.cluster.rds")
#fig0 <- DimPlot(object=atac2, label=F, raster=T)+NoLegend()+
#        theme_bw()+
#        theme(axis.text=element_blank(), 
#              axis.title=element_blank(),
#              axis.ticks=element_blank(), legend.position="none")
#figfn <- "./outs/Figure1.1.1_atac.cluster.png"
#png(figfn, width=200, height=200, res=120)
#print(fig0)
#dev.off()
       
 
}

#####################################
### 2. processing scRNA-seq data  ###
#####################################

if (FALSE){

###
### Read scRNA-seq data
prefix <- "/wsu/home/groups/piquelab/scHOLD/scRNAseq/counts_2021-01-05/"
countfn <- paste(prefix, "scHOLD-RNA-1/outs/filtered_feature_bc_matrix", sep="")
counts <- Read10X(countfn)
sc <- CreateSeuratObject(counts, project="scHOLD")
##
sc <- NormalizeData(sc)
sc <- FindVariableFeatures(sc, nfeatures=3000)
sc <- ScaleData(sc)
sc <- RunPCA(sc, npcs=100)
sc <- RunUMAP(sc, dims=1:50)
sc <- FindNeighbors(sc, dims=1:50, verbose=T)
sc <- FindClusters(sc, resolution=0.15, verbose=T)
###
opfn <- "./outs/1.2_scHOLD.RNAseq.rds"
write_rds(sc, file=opfn)


### summary
### cluster results
fig0 <- DimPlot(object=sc, label=T, raster=F)+NoLegend()+theme_bw()+theme(legend.position="none")
figfn <- "./outs/Figure2.1_cluster.png"
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


####################################
### 3. automatic annotation cell ###
####################################

if (FALSE){

### reference data

ref <- LoadH5Seurat("../pbmc_multimodal.h5seurat")
fig0 <- DimPlot(object=ref, reduction="wnn.umap", group.by="celltype.l1",
                label=T, label.size=3, raster=F, repel=T)+ 
        NoLegend()+
        theme_bw()+
        theme(legend.position="none",
              plot.title=element_text(hjust=0.5))
###
figfn <- "./outs/Figure3.1.ref_umap.png"
png(figfn, width=650, height=500, res=120)
print(fig0)
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



##########################################
### 4. integrate scRNA-seq and sc-ATAC ###
##########################################

if(FALSE){
atac <- read_rds("./outs/1.1_scHOLD.ATAC.cluster.rds")                
sc <- read_rds("./outs/2.2_scHOLD.RNAseq.annot.rds")
sc <- sc%>%
      NormalizeData()%>%
      FindVariableFeatures(nfeatures=3000)%>%
      ScaleData()%>%
      RunPCA(npcs=100)%>%
      RunUMAP(dims=1:50)
       
gene.activities <- GeneActivity(atac, features=VariableFeatures(sc))
atac[["ACTIVITY"]] <- CreateAssayObject(counts=gene.activities)
DefaultAssay(atac) <- "ACTIVITY"
atac <- atac%>%
        NormalizeData()%>%
        ScaleData(features=rownames(sc))        
        
##      
anchors <- FindTransferAnchors(
    reference=sc, query=atac, features=VariableFeatures(object=sc), 
    reference.assay="RNA", query.assay ="ACTIVITY", reduction="cca")
    
pred <- TransferData(anchorset=anchors, refdata=sc$predicted.celltype.l1, 
    weight.reduction=atac[["lsi"]], dims = 2:50)

pred <- pred%>%mutate(barcode=rownames(pred))
x <- atac@meta.data
x <- x%>%left_join(pred,by="barcode")
rownames(x) <- x$barcode
atac <- AddMetaData(atac,metadata=x)
###
opfn <- "./outs/2.1_scHOLD.ATAC.annot.rds"
write_rds(atac, file=opfn)

###
atac <- read_rds("./outs/2.1_scHOLD.ATAC.annot.rds")
fig0 <- DimPlot(atac, reduction="umap.atac", group.by="predicted.id",
              label=T, label.size=3, raster=F, repel=T)+
              NoLegend()+
              theme_bw()+
              theme(legend.position="none",
                    plot.title=element_blank())
figfn <- "./outs/Figure4.1_annot.png"
png(figfn, width=400, height=500, res=120)                                 
print(fig0)
dev.off()

atac <- read_rds("./outs/2.1_scHOLD.ATAC.annot.rds")
fig0 <- DimPlot(atac, reduction="umap.atac", group.by="predicted.id",
              label=T, label.size=3, raster=T, repel=T)+
              NoLegend()+
              theme_bw()+
              theme(legend.position="none", 
              plot.title=element_blank(),
              axis.text=element_blank(), 
              axis.title=element_blank(),
              axis.ticks=element_blank())
figfn <- "./outs/Figure4.1.1_annot.png"
png(figfn, width=400, height=500, res=120)                                 
print(fig0)
dev.off()
}

###
### co-embedding scRNA-seq and scATAC-seq datasets 
if (FALSE){
atac <- read_rds("./outs/2.1_scHOLD.ATAC.annot.rds")
###
sc <- read_rds("./outs/2.2_scHOLD.RNAseq.annot.rds")
sc <- sc%>%
      NormalizeData()%>%
      FindVariableFeatures(nfeatures=3000)%>%
      ScaleData()%>%
      RunPCA(npcs=100)%>%
      RunUMAP(dims=1:50)

### imputation            
genes.use <- VariableFeatures(sc)
refdata <- GetAssayData(sc, assay="RNA", slot="data")[genes.use, ]

anchors <- FindTransferAnchors(
    reference=sc, query=atac, features=VariableFeatures(object=sc), 
    reference.assay="RNA", query.assay ="ACTIVITY", reduction="cca")
     
imputation <- TransferData(anchorset=anchors, refdata=refdata, weight.reduction=atac[["lsi"]], 
    dims = 2:50)
atac[["RNA"]] <- imputation

coembed <- merge(x=sc, y=atac)

# Finally, we run PCA and UMAP on this combined object, to visualize the co-embedding of both
# datasets
coembed <- coembed%>%
           ScaleData(features=genes.use, do.scale=F)%>%
           RunPCA(npcs=100, features=genes.use, verbose=F)%>%
           RunUMAP(dims=1:50)
           
###
opfn <- "./outs/2.3_scHOLD.coembed.rds"
write_rds(coembed, file=opfn)

###
x <- coembed@meta.data
xlabel <- c("scHOLD"="ATAC","SeuratProject"="scRNA-seq")
x$orig.ident2 <- xlabel[x$orig.ident]
coembed@meta.data <- x
fig0 <- DimPlot(coembed, 
        group.by = c("orig.ident2", "predicted.id"),
        reduction="umap", label=T, label.size=3, raster=F, repel=T)+
        NoLegend()&
        theme_bw()+
        theme(legend.position="none",
              plot.title=element_blank())
figfn <- "./outs/Figure5.1_comembed.png"
png(figfn, width=800, height=400, res=120)                                 
print(fig0)
dev.off()     
} ###


###
### directly annotation from CITEseq data

if(FALSE){

ref <- LoadH5Seurat("../pbmc_multimodal.h5seurat")
###
atac <- read_rds("./outs/2.1_scHOLD.ATAC.annot.rds")
DefaultAssay(atac) <- "ACTIVITY"
atac <- SCTransform(atac, assay="ACTIVITY", verbose=FALSE)
anchors <- FindTransferAnchors(reference=ref, query=atac,
           normalization.method="SCT",
           reference.reduction="spca",
           dims=1:50)
###
atac2 <- MapQuery(anchorset=anchors, query=atac, reference=ref,
      refdata = list(
      celltype.l1="celltype.l1",
      celltype.l2="celltype.l2",
      predicted_ADT="ADT"),
      reference.reduction="spca", 
      reduction.model="wnn.umap")
      
opfn <- "./outs/2.1_scHOLD.ATAC.annot2CITEseq.rds"
write_rds(atac2, file=opfn)

### summary results
p1 <- DimPlot(atac2, reduction="ref.umap", group.by="predicted.celltype.l1",
              label=T, label.size=3, raster=F, repel=T)+
              NoLegend()+
              theme_bw()+
              theme(legend.position="none",
                    plot.title=element_text(hjust=0.5))
             
p2 <- DimPlot(atac2, reduction="ref.umap", group.by="predicted.celltype.l2",
              label=T, label.size=3, raster=F, repel=T)+
              NoLegend()+
              theme_bw()+
              theme(legend.position="none",
                    plot.title=element_text(hjust=0.5))
figfn <- "./outs/Figure6.1_ATAC.CITEseq.png"
png(figfn, width=800, height=400, res=120)                                 
print(plot_grid(p1, p2, ncol=2))
dev.off()
###     
}           
           

               