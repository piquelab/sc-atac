#module load R/test_4.0.3
###
###
source("../LibraryPackage.R")
rm(list=ls())

###########################
### Clustering analysis ###
###########################

if(FALSE){
atac <- read_rds("./outs/1_seurat.merge.rds")
                
atac2 <- RunTFIDF(atac)
atac2 <- FindTopFeatures(atac2, min.cutoff="q0")
atac2 <- RunSVD(atac2,n=100)

fig1 <- DepthCor(atac2)+
        theme(plot.title=element_text(size=10),plot.subtitle=element_text(size=10))
figfn <- "./outs/Figure4.depthcor.png"
png(figfn, width=600, height=400, res=120)
print(fig1)
dev.off()

##
atac2 <- RunUMAP(atac2, reduction="lsi",dims=2:50)
atac2 <- FindNeighbors(atac2, reduction="lsi", dims=2:50)
atac2 <- FindClusters(atac2, resolution=0.15, verbose=FALSE, algorithm=3)

write_rds(atac2, file="./outs/2_seurat.cluster.rds")
}

###############
### summary ###
###############

if(TRUE){
atac2 <- read_rds("./outs/2_seurat.cluster.rds")

fig2 <- DimPlot(object=atac2, label=TRUE, raster=F)+NoLegend()+theme_bw()+theme(legend.position="none")
figfn <- "./outs/Figure5.1_cluster.png"
png(figfn, width=500, height=400, res=120)
print(fig2)
dev.off()

atac2$Batch <- gsub("-.*", "", colnames(atac2))
fig3 <- DimPlot(object=atac2, group.by="Batch", raster=F)+NoLegend()+
        facet_grid(~Batch)+ 
        theme_bw()+
        theme(legend.position="none",
              plot.title=element_blank())
figfn <- "./outs/Figure5.2_BATCH.png"
png(figfn, width=700, height=350, res=120)
print(fig3)
dev.off()

col1 <- c("e"="#828282", 
           "a"="#fb9a99", "b"="#e31a1c",
           "c"="#a6cee3", "d"="#1f78b4")
atac2$treat <- gsub(".*-ATAC-|_.*", "", colnames(atac2))
alpha <- c("LPS"="a", "LPS-DEX"="b", "PHA"="c", "PHA-DEX"="d", "CTRL"="e")
atac2$treat.alpha <- alpha[atac2$treat]

fig3 <- DimPlot(object=atac2, group.by="treat.alpha", raster=F)+NoLegend()+
        facet_wrap(~treat.alpha, ncol=2, 
                   labeller=as_labeller(c("a"="LPS", "b"="LPS-DEX", "c"="PHA", "d"="PHA-DEX", "e"="CTRL")))+
        scale_colour_manual(values=col1)+
        theme_bw()+
        theme(legend.position="none",
              plot.title=element_blank())
figfn <- "./outs/Figure5.3_treat.png"
png(figfn, width=650, height=650, res=120)
print(fig3)
dev.off()

}


                