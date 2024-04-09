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
library(cicero, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
library(monocle3)
library(EnsDb.Hsapiens.v75)
###
library(ggplot2)
library(patchwork)
library(cowplot)
library(RColorBrewer)
library(viridis)
library(ggrastr)
theme_set(theme_grey())

rm(list=ls())

###
###
outdir <- "./6_pub.outs/1_main_plots/"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)

## outdir <- "./6_pub.outs/poster/"
## if (!file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)




############
### UMAP ###
############

fn <- "./4.2_Integrate.outs/3_scATAC.annot.rds"
atac <- read_rds(fn)
## DefaultAssay(atac) <- "ACTIVITY"

x <- atac@meta.data
MCl_lab <- c("Bcell"="B cell", "Monocyte"="Monocyte", "NKcell"="NK cell", "Tcell"="T cell", "DC"="DC")
x <- x%>%mutate(MCl2=MCl_lab[as.character(MCls)])

atac <- AddMetaData(atac, metadata=x)

col2 <- c("#4daf4a", "#828282", "#984ea3", "#aa4b56", "#ffaa00")

## col2 <- c("Bcell"="#4daf4a", "Monocyte"="#984ea3", "NKcell"="#aa4b56", "Tcell"="#ffaa00", "DC"="#828282")
          
p0 <- DimPlot(atac, reduction="umap.atac", group.by="MCl2", cols=col2, label=T, label.size=5, raster=F)+
   theme_bw()+   
   xlab("UMAP_1")+ylab("UMAP_2")+ 
   theme(##legend.position="none",
         plot.title=element_blank(),
         axis.text=element_text(size=12),
         axis.title=element_text(size=12))

   ## ## ## guides(col=guide_legend(override.aes=list(size=2),ncol=3))+
   ## theme(legend.title=element_blank(),
   ##       legend.key.size=grid::unit(0.8,"lines"))
figfn <- paste(outdir, "Figure1.1_umap_atac.MCls.png", sep="")
png(figfn, width=480, height=500, res=120)
print(p0)
dev.off()




### poster
          
p0 <- DimPlot(atac, reduction="umap.atac", group.by="MCl2", cols=col2, label=F, label.size=5, raster=F)+
   theme_bw()+   
   xlab("UMAP_1")+ylab("UMAP_2")+   
   guides(col=guide_legend(override.aes=list(size=2),ncol=1))+ 
   theme(legend.position="none",
         plot.title=element_blank(),
         axis.text=element_text(size=12),
         axis.title=element_text(size=12),
         legend.title=element_blank(),
         legend.key.size=grid::unit(0.8,"lines"))

figfn <- paste(outdir, "Figure1.1_umap.MCls.png", sep="")
png(figfn, width=420, height=520, res=120)
print(p0)
dev.off()



###########################
#### coverage plots
##########################

rm(list=ls())

outdir <- "./6_pub.outs/1_main_plots/"

source("annot_plot.R")

fn <- "./4.2_Integrate.outs/3_scATAC.annot.rds"
atac <- read_rds(fn)
## DefaultAssay(atac) <- "ACTIVITY"


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
  s0 <- ifelse(i==1, 0.5e+03, 1e+03)  
  region <- FindRegion(atac, region=geneId, assay="ATAC", extend.upstream=s0, extend.downstream=1e+03)
  p2 <- AnnotationPlot(atac, region=region)&
     theme(axis.title.x=element_blank(),
           axis.text.x=element_blank(),
           axis.ticks.x=element_blank(),
           axis.title.y=element_blank())
  ##  
  ll <- length(p2$layers)
   p2$layers[[ll]]$aes_params$size <- 3
  p2$layers[[ll]]$aes_params$fontface <- "italic"  
  ###  
  p <- wrap_plots(p1, p2, ncol=1, heights=c(8, 1.2))  
  ##        
  p
})

###
figfn <- paste(outdir, "Figure1c_MCls.coverage.png", sep="")
png(figfn, width=700, height=360, res=100)
plot_grid(plotlist=figs_ls, ncol=5)
dev.off()
