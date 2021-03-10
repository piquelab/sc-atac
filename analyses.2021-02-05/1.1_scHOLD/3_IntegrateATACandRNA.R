##
source("../LibraryPackage.R")

outdir <- "./3_Integrate.outs/"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)
rm(list=ls())

#######################################
### integrate scRNA-seq and sc-ATAC ###
#######################################

if(TRUE){
               
sc <- read_rds("./outs/2.2_scHOLD.RNAseq.annot.rds")
sc <- sc%>%
      NormalizeData()%>%
      FindVariableFeatures(nfeatures=3000)%>%
      ScaleData()%>%
      RunPCA(npcs=100)%>%
      RunUMAP(dims=1:50)

atac <- read_rds("./outs/1.1_scHOLD.ATAC.cluster.rds")       
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

### Using predicted.celltype.l1     
pred <- TransferData(anchorset=anchors, refdata=sc$predicted.celltype.l1, 
    weight.reduction=atac[["lsi"]], dims = 2:50)

pred <- pred%>%mutate(barcode=rownames(pred))
x <- atac@meta.data
x <- x%>%left_join(pred,by="barcode")
rownames(x) <- x$barcode
atac2 <- AddMetaData(atac,metadata=x)
###
opfn <- "./3_Integrate.outs/1.1_scHOLD.ATAC.annot1.rds"
write_rds(atac2, file=opfn)    

###
### Using predicted.celltype.l2    
pred <- TransferData(anchorset=anchors, refdata=sc$predicted.celltype.l2, 
        weight.reduction=atac[["lsi"]], dims = 2:50)

pred <- pred%>%mutate(barcode=rownames(pred))
x <- atac@meta.data
x <- x%>%left_join(pred,by="barcode")
rownames(x) <- x$barcode
atac2 <- AddMetaData(atac,metadata=x)
###
opfn <- "./3_Integrate.outs/1.2_scHOLD.ATAC.annot2.rds"
write_rds(atac2, file=opfn)
}###


###Summary
if (FALSE){
###
atac2 <- read_rds("./3_Integrate.outs/1.1_scHOLD.ATAC.annot1.rds")
p1 <- DimPlot(object=atac2, reduction="umap.atac", label=T, raster=F)+NoLegend()+
      theme_bw()+
      theme(legend.position="none")
p2 <- DimPlot(atac2, reduction="umap.atac", group.by="predicted.id",
              label=T, label.size=3, raster=F, repel=T)+
              NoLegend()+
              theme_bw()+
              theme(plot.title=element_blank(),
                    legend.key.size=grid::unit(1,"lines"))
figfn <- "./3_Integrate.outs/Figure1.1_annot.png"
png(figfn, width=800, height=400, res=120)                                 
print(plot_grid(p1, p2, ncol=2, rel_widths=c(1,1.3)))
dev.off()

atac2 <- read_rds("./3_Integrate.outs/1.2_scHOLD.ATAC.annot2.rds")
p3 <- DimPlot(atac2, reduction="umap.atac", group.by="predicted.id",
              label=T, label.size=2.5, raster=F, repel=T)+
              NoLegend()+
              theme_bw()+
              theme(plot.title=element_blank(),
                    legend.text=element_text(size=8),
                    legend.key.size=grid::unit(0.5,"lines"))
figfn <- "./3_Integrate.outs/Figure1.2_annot.png"
png(figfn, width=650, height=400, res=120)                                 
print(p3)
dev.off()


#atac <- read_rds("./outs/2.1_scHOLD.ATAC.annot.rds")
#fig0 <- DimPlot(atac, reduction="umap.atac", group.by="predicted.id",
#              label=T, label.size=3, raster=T, repel=T)+
#              NoLegend()+
#              theme_bw()+
#              theme(legend.position="none", 
#              plot.title=element_blank(),
#              axis.text=element_blank(), 
#              axis.title=element_blank(),
#              axis.ticks=element_blank())
#figfn <- "./outs/Figure4.1.1_annot.png"
#png(figfn, width=400, height=500, res=120)                                 
#print(fig0)
#dev.off()

### features
#x0 <- c("MS4A1", "CD79A", "MS4A7", "CD14", "GNLY", "NKG7", "CD3D", "CD8A")
#atac2 <- read_rds("./outs/2.1_scHOLD.ATAC.annot.rds")
#fig0 <- FeaturePlot(object=atac2, features=x0, ncol=4, raster=F)&
#        theme_bw()+
#           theme(legend.title=element_blank(),
#                 legend.key.size=grid::unit(0.5,"lines"),
#                 plot.title=element_text(size=12, hjust=0.5))
#png("./3_Integrate.outs/Figure2.1_feature.png", width=950, height=500, res=100)
#print(fig0)
#dev.off()
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
### UMAP from atac
fn <- "./outs/2.1_scHOLD.ATAC.annot2CITEseq.rds"
atac2 <- read_rds(fn)
p1 <- DimPlot(atac2, reduction="umap.atac", group.by="predicted.celltype.l1",
              label=T, label.size=3, raster=F, repel=T)+
              theme_bw()+
              theme(plot.title=element_blank())
             
figfn <- "./outs/Figure6.2_ATAC.CITEseq.atacUMAP.png"
png(figfn, width=500, height=500, res=120)                                 
print(p1)
dev.off()

}           
           

               