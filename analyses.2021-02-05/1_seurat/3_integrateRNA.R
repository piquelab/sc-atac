#module load R/test_4.0.3
###
###
source("../LibraryPackage.R")
rm(list=ls())

###
atac <- read_rds("./outs/2_seurat.cluster.rds")
gene.activities <- GeneActivity(atac)

atac[["RNA"]] <- CreateAssayObject(counts=gene.activities)
atac2 <- NormalizeData(object=atac,
                       assay="RNA",
                       normalization.method="LogNormalize")

                       
####################                       
## summary plots ###
####################

DefaultAssay(atac2) <- "RNA"

x0 <- c("MS4A1", "CD79A", "MS4A7", "CD14", "GNLY", "NKG7", "CD3D", "CD8A")

fig1 <- FeaturePlot(object=atac2, 
                    features=x0, ncol=4, raster=F)&
        #scale_color_gradient("",low="lightgrey",high="blue")+
        theme_bw()+
           theme(legend.title=element_blank(),
                 legend.key.size=grid::unit(0.5,"lines"),
                 plot.title=element_text(size=12, hjust=0.5))
png("./outs/Figure6.1_MCls.feature.png", width=950, height=500, res=100)
print(fig1)
dev.off()   
                     
                        