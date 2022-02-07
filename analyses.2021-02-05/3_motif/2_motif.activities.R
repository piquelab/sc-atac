##
library(tidyverse)
library(Seurat)
library(parallel)
library(SeuratDisk)
library(SeuratData)
library(SeuratObject)
library(Signac)
library(SeuratWrappers)
library(SeuratData)
library(chromVAR)
    
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Hsapiens.1000genomes.hs37d5, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
library(ChIPseeker, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
###
library(ggplot2)
library(cowplot, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
library(grid)
library(gridExtra)
library(ggExtra)
library(RColorBrewer)
library(ggsci)
library(viridis)
library(ComplexHeatmap)
library(circlize)
theme_set(theme_grey())

###
outdir <- "./2_motif.activities.outs/"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)

##
### calculate motif
fn <- "./1_motif.outs/1_scATAC.motif.rds"
atac <- read_rds(fn)

atac <- RunChromVAR(object=atac,
   genome=BSgenome.Hsapiens.1000genomes.hs37d5)

opfn <- "./2_motif.activities.outs/1_scATAC.motifActivities.rds"
write_rds(atac, opfn)


### tsne plots
## atac <- read_rds("./2_motif.activities.outs/1_scATAC.motifActivities.rds")
## X <- atac@assays$chromvar@data

## library(tsne)
## xtsne <- tsne(t(X))
## write_rds(xtsne, "./2_motif.activities.outs/1.1_tsne.rds")


###
### differential motif activities
atac <- read_rds("./2_motif.activities.outs/1_scATAC.motifActivities.rds")
atac2 <- subset(atac, subset=MCls!="DC")

Y <- atac2@assays$chromvar@data
meta <- atac2@meta.data%>%
   mutate(treat=gsub(".*-ATAC-|_.*", "", NEW_BARCODE),
          bti=paste(MCls, treat, SNG.BEST.GUESS, sep="_"))%>%
   dplyr::select(NEW_BARCODE, bti) 

###
dd <- meta%>%
   group_by(bti)%>%
   summarise(ncell=n(), .groups="drop")%>%
   ungroup()%>%as.data.frame()

###
bti <- factor(meta$bti)
X <- model.matrix(~0+bti)
YtX <- Y %*% X
YtX <- as.matrix(YtX)
colnames(YtX) <- gsub("^bti", "", colnames(YtX))

YtX_ave <- sweep(YtX, 2, dd$ncell, "/")

dd0 <- dd%>%dplyr::filter(ncell>20)

YtX_sel <- YtX_ave[,dd0$bti]

write_rds(YtX_sel, "./2_motif.activities.outs/2_motif.ave.rds")


###
### heatmap
mat <- read_rds("./2_motif.activities.outs/2_motif.ave.rds")
meta <- str_split(colnames(mat), "_", simplify=T)%>%as.data.frame()
names(meta) <- c("MCls", "treat", "sampleID")
###
b <- as.vector(mat)
breaks <- quantile(b, probs=seq(0, 1, length.out=100))
col_fun <-  colorRamp2(breaks,
   colorRampPalette(rev(brewer.pal(n=7, name="RdBu")))(100))
column_ha <- HeatmapAnnotation(
   celltype=meta$MCls,
   treatment=meta$treat,
   col=list(
      celltype=c("Bcell"="#4daf4a", "Monocyte"="#984ea3",
                 "NKcell"="#aa4b56", "Tcell"="#ffaa00"),
      treatment=c("CTRL"="#828282", "LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
                  "PHA"="#a6cee3", "PHA-DEX"="#1f78b4")))

fig <- Heatmap(mat, col=col_fun,
   cluster_rows=T, cluster_columns=T,
   show_row_dend=T, show_column_dend=T,
   top_annotation=column_ha,
   heatmap_legend_param=list(title="motif activities",
      title_gp=gpar(fontsize=10),
      labels_gp=gpar(fontsize=10)),
   show_row_names=F, show_column_names=F,
   use_raster=F, raster_device="png")

figfn <- "./2_motif.activities.outs/Figure1.1_heatmap.motif.activities.png"
png(figfn, height=800, width=900, res=120)
set.seed(0)
fig <- draw(fig)
dev.off()


###############################
### Differential activities ###
###############################
### myDE
myDE <- function(y, X, gene){

   con.ls <- list("LPS"=c("CTRL","LPS"),
      "LPS-DEX"=c("LPS","LPS-DEX"),
      "PHA"=c("CTRL", "PHA"),
      "PHA-DEX"=c("PHA","PHA-DEX"))
   Contrast <- c("LPS", "LPS-DEX", "PHA", "PHA-DEX")

   ## linear regression
   x1 <- X[,1]
   lm0 <- try(lm(y~0+x1), silent=T)

   if ( (class(lm0)!="try-error")){
      b <- coef(lm0)
      nn <- gsub("x[12]", "", names(b))
      names(b) <- nn
      vb <- diag(vcov(lm0))
      names(vb) <- nn

      ## Contrast
      dd <- lapply(Contrast,function(one){
        con0 <- con.ls[[one]]
        if ( all(con0%in%nn)){
          bhat <- b[con0[2]]-b[con0[1]]
          sdhat <- sqrt(vb[con0[2]]+vb[con0[1]])
          z <- bhat/sdhat
          p <- 2*pnorm(-abs(z))
          dd1 <- data.frame(gene=gene, beta=bhat, stderr=sdhat, pval=p,
                            contrast=one)
        }else{
          dd1 <- NA                                                                       }
        dd1
      })
      dd <- dd[!is.na(dd)]
      dd <- do.call(rbind,dd)
      ###
   }else{
      dd <- NA
   } ## End if try-error
   dd 
} ###



mat <- read_rds("./2_motif.activities.outs/2_motif.ave.rds")
bti <-colnames(mat)
cvt0 <- str_split(bti, "_", simplify=T)%>%as.data.frame()
cvt <- data.frame(bti=bti, MCls=cvt0[,1], treat=cvt0[,2], sampleID=cvt0[,3])

### Differential procedure
MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
res <- map_dfr(MCls, function(oneMCl){

   cvti <- cvt%>%dplyr::filter(MCls==oneMCl)
   mati <- mat[,cvti$bti]
   X <- data.frame(x1=cvti$treat)
   rn <- rownames(mati)
   TMP <- mclapply(rn, function(ii){
     y <- mati[ii,]
     dd <- myDE(y, X, ii)
     dd
   },mc.cores=1)
   TMP <- TMP[!is.na(TMP)]
   TMP <- do.call(rbind, TMP)%>%as.data.frame()%>%mutate(MCls=oneMCl)
   TMP 
})    

### add qvalue
res2 <- res%>%group_by(MCls, contrast)%>%
   mutate(qval=p.adjust(pval, "BH"))%>%
   ungroup()%>%
   as.data.frame()

atac <- read_rds("./2_motif.activities.outs/1_scATAC.motifActivities.rds")
motif <- Motifs(atac)
x <- unlist(motif@motif.names)

res2 <- res2%>%mutate(motif=x[gene])
    
###
opfn <- "./2_motif.activities.outs/3_motif.diff.results.rds"
write_rds(res2, opfn)


###
### significant motif
res <- read_rds("./2_motif.activities.outs/3_motif.diff.results.rds")
sigs <- res%>%dplyr::filter(qval<0.1,abs(beta)>1.41)%>%
   group_by(MCls, contrast)%>%
   summarise(ny=n(), .groups="drop")

topmotif <- res%>%
   dplyr::filter(qval<0.1, abs(beta)>1.41)%>%
   dplyr::pull(motif)
topmotif <- sort(unique(topmotif))

### data for heatmap
res <- res%>%mutate(condition=paste(MCls, contrast, sep="_"))
condition <- sort(unique(res$condition))
mat <- map_dfc(condition, function(ii){
   res2 <- res%>%dplyr::filter(condition==ii)
   b <- res2$beta
   names(b) <- res2$motif
   b[topmotif]
})
mat <- as.matrix(mat)
colnames(mat) <- condition
rownames(mat) <- topmotif

###
###
b <- as.vector(mat)
breaks <- quantile(b, probs=seq(0, 1, length.out=100))
col_fun <-  colorRamp2(breaks,
   colorRampPalette(rev(brewer.pal(n=7, name="RdBu")))(100))
column_ha <- HeatmapAnnotation(
   celltype=gsub("_.*", "", condition),
   treatment=gsub(".*_", "", condition),
   col=list(
      celltype=c("Bcell"="#4daf4a", "Monocyte"="#984ea3",
                 "NKcell"="#aa4b56", "Tcell"="#ffaa00"),
      treatment=c("LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
                  "PHA"="#a6cee3", "PHA-DEX"="#1f78b4")))

fig <- Heatmap(mat, col=col_fun,
   cluster_rows=T, cluster_columns=T,
   show_row_dend=F, show_column_dend=T,
   top_annotation=column_ha,
   heatmap_legend_param=list(title="Diff motif",
      title_gp=gpar(fontsize=10),
      labels_gp=gpar(fontsize=10)),
   show_row_names=T, show_column_names=T,
   row_names_side="left",
   column_names_gp=gpar(fontsize=9),
   row_names_gp=gpar(fontsize=6),
   use_raster=F, raster_device="png")

figfn <- "./2_motif.activities.outs/Figure1.2_heatmap.motif.activities.png"
png(figfn, height=800, width=700, res=120)
set.seed(0)
fig <- draw(fig)
dev.off()


##
#################################
### write differential motif ###
################################

res <- read_rds("./2_motif.activities.outs/3_motif.diff.results.rds")
##
## motif.name <- sort(unique(res$motif))
## motif.id <- ConvertMotifID(object=motif, name=motif.name)
## motif.annot <- data.frame(motif.name=motif.name, motif.id=motif.id)

## res <- res%>%left_join(motif.annot, by=c("motif"="motif.name"))

sigs <- res%>%mutate(condition=paste(MCls, contrast, sep="_"))%>%
   dplyr::filter(qval<0.1, beta>0)%>%  ### abs(beta)>0.32)%>% ##90%
   group_by(condition)%>%
   summarise(ny=n(), .groups="drop")

###
###
res2 <- res%>%mutate(condition=paste(MCls, contrast, sep="_"))%>%
   dplyr::filter(qval<0.1, beta>0)##abs(beta)>0.32)

condition <- sort(unique(res2$condition))
motif_response <- map(condition, function(one){
    res2%>%dplyr::filter(condition==one)%>%dplyr::pull(gene)
})
names(motif_response) <- condition

opfn <- paste(outdir, "4.2_response.motif.positive.rds", sep="")
write_rds(motif_response, opfn)


#####################
### resonse motif ###
#####################

outdir <- "./2_motif.activities.outs/"

res <- read_rds("./2_motif.activities.outs/3_motif.diff.results.rds")

##
oneMCl <- "Tcell"
res2 <- res%>%filter(MCls==oneMCl)

th0 <- quantile(abs(res2$beta),probs=0.9)
th0 <- 0.21
### CTRL
res3 <- res2%>%filter(contrast=="LPS"|contrast=="PHA", qval<0.1, abs(beta)>th0)%>%
    arrange(desc(abs(beta)))
opfn <- paste(outdir, "5_",  oneMCl, "_CTRL.motif.txt", sep="")
write.table(unique(res3$gene), opfn, row.names=F, quote=F, col.names=F)


### LPS-EtOH
res3 <- res2%>%filter(contrast=="LPS", qval<0.1, abs(beta)>th0)%>%
    arrange(desc(abs(beta)))
opfn <- paste(outdir, "5_",  oneMCl, "_LPS-EtOH.motif.txt", sep="")
write.table(unique(res3$gene), opfn, row.names=F, quote=F, col.names=F)


### LPS-DEX
res3 <- res2%>%filter(contrast=="LPS-DEX", qval<0.1, abs(beta)>th0)%>%
    arrange(desc(abs(beta)))
opfn <- paste(outdir, "5_",  oneMCl, "_LPS-DEX.motif.txt", sep="")
write.table(unique(res3$gene), opfn, row.names=F, quote=F, col.names=F)
 

### PHA-EtOH
res3 <- res2%>%filter(contrast=="PHA", qval<0.1, abs(beta)>th0)%>%
    arrange(desc(abs(beta)))
opfn <- paste(outdir, "5_",  oneMCl, "_PHA-EtOH.motif.txt", sep="")
write.table(unique(res3$gene), opfn, row.names=F, quote=F, col.names=F)


### PHA-DEX
res3 <- res2%>%filter(contrast=="PHA-DEX", qval<0.1, abs(beta)>th0)%>%
    arrange(desc(abs(beta)))
opfn <- paste(outdir, "5_",  oneMCl, "_PHA-DEX.motif.txt", sep="")
write.table(unique(res3$gene), opfn, row.names=F, quote=F, col.names=F)







##############################################################
### identify motif with high activities in some celll type ###
##############################################################
atac2 <- subset(atac2, subset=MCls!="DC")

Y <- atac2@assays$chromvar@data
meta <- atac2@meta.data%>%
   mutate(treat=gsub(".*-ATAC-|_.*", "", NEW_BARCODE),
          bti=paste(MCls, treat, SNG.BEST.GUESS, sep="_"))%>%
   dplyr::select(NEW_BARCODE, bti) 
