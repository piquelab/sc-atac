#module load R/test_4.0.3
###
###
source("../LibraryPackage.R")
rm(list=ls())

atac <- read_rds("./outs/1_seurat.merge.rds")
atac2 <- subset(x=atac,
                subset=peak_region_fragments>3000&
                peak_region_fragments<25000&
                pct_reads_in_peaks>20&
                blacklist_ratio<0.2&
                nucleosome_signal<10&
                TSS.enrichment>2)
                
atac2 <- RunTFIDF(atac2)
atac2 <- FindTopFeatures(atac2, min.cutoff="q0")
atac2 <- RunSVD(atac2)

fig1 <- DepthCor(atac2)
figfn <- "./outs/Figure4.depthcor.png"
png(figfn, width=600, height=400, res=120)
print(fig1)
dev.off()

##
atac2 <- RunUMAP(atac2,,reduction="lsi",dims=2:30)
atac2 <- FindNeighbors(atac2, reduction="lsi", dims=2:30)
atac2 <- FindClusters(atac2, verbose=FALSE, algorithm=3)

write_rds(atac2, file="./outs/2_seurat.cluster.rds")

###
###summary
atac2 <- read_rds("./outs/2_seurat.cluster.rds")

fig2 <- DimPlot(object=atac2, label=TRUE)+NoLegend()+theme_bw()+theme(legend.position="none")
figfn <- "./outs/Figure5.1_cluster.png"
png(figfn, width=600, height=400, res=120)
print(fig2)
dev.off()


atac2$EXP <- gsub("_.*", "", colnames(atac2))
fig3 <- DimPlot(object=atac2, group.by="EXP")+NoLegend()+theme_bw()+theme(plot.title=element_blank())
figfn <- "./outs/Figure5.2_EXP.png"
png(figfn, width=800, height=500, res=120)
print(fig3)
dev.off()


                