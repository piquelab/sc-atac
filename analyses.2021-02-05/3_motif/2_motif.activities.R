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

rm(list=ls())

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
   x2 <- X[,2] 
   lm0 <- try(lm(y~0+x1+x2), silent=T)

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
   X <- data.frame(x1=cvti$treat, x2=cvti$sampleID)
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
motif.name <- x
motif.id <- names(x)
names(motif.id) <- motif.name

res2 <- res2%>%mutate(motif=motif.id[gene])
    
###
opfn <- "./2_motif.activities.outs/3.2_indi_motif.diff.results.rds"
write_rds(res2, opfn)




###
### Heatmap for publication 



########################################
### show heatmap ofsignificant motif ###
########################################

res <- read_rds("./2_motif.activities.outs/3_motif.diff.results.rds")


###
## fn <- "./1.2_motif.outs/3_motif.enrich.direction.rds"
## enrich <- read_rds(fn)
## drt2 <- c("0"="Down", "1"="Up")
## enrich <- enrich%>%mutate(direction2=drt2[as.character(direction)],
##    cluster=paste(contrast, MCls, direction2, sep="."),
##    newCluster=clst2[cluster])
###
## topmotif <- enrich%>%
##    group_by(cluster)%>%
##    top_n(n=6, wt=fold.enrichment)%>%ungroup()%>%
##     dplyr::pull(motif)%>%unique()

## topmotif2 <- enrich%>%
##     filter(motif%in%topmotif, fold.enrichment>1.41, qvalue.fisher<0.1)%>%
##     dplyr::pull(motif.name)%>%unique()

sigs <- res%>%dplyr::filter(qval<0.1)%>% ##,abs(beta)>1.41)%>%
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
   b <- res2$beta ##/res2$stderr
   names(b) <- res2$motif
   b[topmotif]
})
mat <- as.matrix(mat)
colnames(mat) <- condition
rownames(mat) <- topmotif
 
###
###
res2 <- res%>%dplyr::filter(motif%in%topmotif)%>%mutate(is_sig=ifelse(qval<0.1, 1, 0))
imat <- res2%>%pivot_wider(id_cols=motif, names_from=condition, values_from=is_sig)
imat2 <- as.matrix(imat[,-1])
rownames(imat2) <- as.character(imat$motif)
imat2 <- imat2[topmotif,]



## b <- as.vector(mat)
mat2 <- mat*imat2
colnames(mat2) <- gsub("-", "+", colnames(mat2))


condition2 <- c("Bcell_LPS", "Bcell_PHA", "Monocyte_LPS", "Monocyte_PHA",
                "NKcell_LPS", "NKcell_PHA", "Tcell_LPS", "Tcell_PHA",
               "Bcell_LPS+DEX", "Bcell_PHA+DEX", "Monocyte_LPS+DEX", "Monocyte_PHA+DEX",
                "NKcell_LPS+DEX", "NKcell_PHA+DEX", "Tcell_LPS+DEX", "Tcell_PHA+DEX")

## condition2 <- c("Bcell_LPS", "Bcell_PHA", "Bcell_LPS+DEX", "Bcell_PHA+DEX",
##                  "Monocyte_LPS", "Monocyte_PHA","Monocyte_LPS+DEX", "Monocyte_PHA+DEX",
##                  "NKcell_LPS", "NKcell_PHA", "NKcell_LPS+DEX", "NKcell_PHA+DEX",
##                  "Tcell_LPS", "Tcell_PHA", "Tcell_LPS+DEX", "Tcell_PHA+DEX")

## condition2 <- c("Bcell_LPS", "Monocyte_LPS", "NKcell_LPS", "Tcell_LPS",  
##                 "Bcell_PHA", "Monocyte_PHA", "NKcell_PHA", "Tcell_PHA",  
##                 "Bcell_LPS+DEX", "Monocyte_LPS+DEX", "NKcell_LPS+DEX", "Tcell_LPS+DEX",
##                 "Bcell_PHA+DEX", "Monocyte_PHA+DEX", "NKcell_PHA+DEX", "Tcell_PHA+DEX")

mat2 <- mat2[, condition2]


b <- as.vector(mat2)
b0 <- quantile(b[b<0], probs=seq(0, 1, length.out=49))
b1 <- 0
b2 <- quantile(b[b>0], probs=seq(0, 1, length.out=49))
breaks <- c(b0, b1, b2)

## breaks <- quantile(b, probs=seq(0, 1, length.out=100),na.rm=T)
 
col_fun <-  colorRamp2(breaks,
   colorRampPalette(rev(brewer.pal(n=7, name="RdBu")))(99))

column_ha <- HeatmapAnnotation(
   celltype=gsub("_.*", "", condition2),
   contrast=gsub("-", "+", gsub(".*_", "", condition2)),
   col=list(
      celltype=c("Bcell"="#4daf4a", "Monocyte"="#984ea3",
                 "NKcell"="#aa4b56", "Tcell"="#ffaa00"),
      contrast=c("LPS"="#fb9a99", "LPS+DEX"="#e31a1c",
                  "PHA"="#a6cee3", "PHA+DEX"="#1f78b4")),
   annotation_legend_param=list(celltype=list(labels_gp=gpar(fontsize=8)),
                                contrast=list(labels_gp=gpar(fontsize=8))))


###
### row annotation
 
### cluster pattern
fn <- "./2_motif.activities.outs/Figure1.4_row_cluster.txt"
geneCL <- read.table(fn, header=T)
geneCL <- geneCL%>%mutate(cluster=paste("cluster", cluster, sep=""))

df <- data.frame(cluster=c("cluster1", "cluster2", "cluster3", "cluster4"),
                 cluster2=c("4", "1", "2", "3"))
geneCL <- geneCL%>%left_join(df, by="cluster")

 
anno_df <- data.frame(motif=rownames(mat2))
anno_df2 <- anno_df%>%left_join(geneCL[,2:3], by="motif")%>%dplyr::select(Pattern=cluster2)


row_ha <- rowAnnotation(df=anno_df2,
   col=list(Pattern=c("1"="#1b9e77", "2"="#d95f02",
                     "3"="#e7298a", "4"="#7570b3")),
   annotation_legend_param=list(
   Pattern=list(labels_gp=gpar(fontsize=8), title_gp=gpar(fontsize=10), grid_width=grid::unit(0.6, "cm"),
                grid_height=grid::unit(0.8, "cm"))),
   show_annotation_name=F, simple_anno_size=grid::unit(0.5, "cm")) 

                        

fig <- Heatmap(mat2, col=col_fun,
   cluster_rows=T, cluster_columns=F,
   show_row_dend=T,  show_column_dend=F,
   top_annotation=column_ha,
   right_annotation=row_ha,
   heatmap_legend_param=list(title="Diff motif",
      title_gp=gpar(fontsize=10, font=2),
      labels_gp=gpar(fontsize=10),
      legend_height=grid::unit(4, "cm"),
      grid_width=grid::unit(0.4, "cm")),
   show_row_names=T, show_column_names=T,
   ##row_names_side="left",
   column_names_gp=gpar(fontsize=10),
   row_names_gp=gpar(fontsize=8),
   use_raster=F, raster_device="png")
  
figfn <- "./2_motif.activities.outs/Figure1.4_heatmap.motif.activities.png"
png(figfn, height=850, width=900, res=120)
set.seed(0)
fig <- draw(fig)
dev.off()


###
### cluster outputs
hmap <- Heatmap(mat2, cluster_rows=T, cluster_columns=F, row_split=4)
set.seed(0)
hmap <- draw(hmap)
cl <- row_order(hmap)


DF_clu <- NULL
for (i in 1:length(cl)){
   ## 
   tmp <- data.frame(cluster=i, motif=topmotif[cl[[i]]])
   DF_clu <- rbind(DF_clu,tmp)
}

opfn <- paste(outdir, "Figure1.4_row_cluster.txt", sep="")
write.table(DF_clu, file=opfn, row.names=F, col.names=T, quote=F, sep="\t")


### swap
## row_ha <- HeatmapAnnotation(
##    celltype=gsub("_.*", "", condition),
##    treatment=gsub(".*_", "", condition),
##    col=list(
##       celltype=c("Bcell"="#4daf4a", "Monocyte"="#984ea3",
##                  "NKcell"="#aa4b56", "Tcell"="#ffaa00"),
##       treatment=c("LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
##                   "PHA"="#a6cee3", "PHA-DEX"="#1f78b4")))

## fig <- Heatmap(t(mat), col=col_fun,
##    cluster_rows=T, cluster_columns=T,
##    show_row_dend=T, show_column_dend=F,
##    top_annotation=row_ha,
##    heatmap_legend_param=list(title="Diff motif",
##       title_gp=gpar(fontsize=10),
##       labels_gp=gpar(fontsize=10), legend_direction="horizontal"),
##    show_row_names=T, show_column_names=T,
##    row_names_side="left",
##    column_names_gp=gpar(fontsize=6),
##    row_names_gp=gpar(fontsize=9),
##    use_raster=F, raster_device="png")

## figfn <- "./2_motif.activities.outs/Figure1.3_heatmap.motif.activities.png"
## png(figfn, height=550, width=700, res=120)
## set.seed(0)
## fig <- draw(fig)
## dev.off()


###
### End


### output motifs




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




##################################
### define  response motif (1) ###
##################################

rm(list=ls())

outdir <- "./2_motif.activities.outs/MotifList/"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)


res <- read_rds("./2_motif.activities.outs/3_motif.diff.results.rds")
res <- res%>%mutate(comb=paste(MCls, contrast, sep="_"))

###########################################
### conditions separately define motifs ###
###########################################

conditions <- sort(unique(res$comb))
for (ii in conditions){
   ###
   oneMCl <- gsub("_.*", "", ii)
   res2 <- res%>%filter(MCls==oneMCl)
    
   if( oneMCl!="Monocyte"){ 
      th0 <- quantile(abs(res2$beta), probs=0.9)
   }else{   
      th0 <- quantile(abs(res$beta),probs=0.9)
   }
    
   ##
   motifs <- res%>%filter(comb==ii, qval<0.1, abs(beta)>th0)%>%pull(gene)%>%unique()
   opfn <- paste(outdir, "1_comb_", ii, "_motif.txt", sep="") 
   write.table(motifs, opfn, row.names=F, quote=F, col.names=F)
   ## 
   cat(ii, length(motifs), "\n")
}



#######################################################
### for contrast union of motifs across cell-types  ###
#######################################################
 
contrast <- sort(unique(res$contrast))
MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
for (ii in contrast){
    
   ###
   motifs <- lapply(MCls, function(oneMCl){
   ##
     comb <- paste(oneMCl, ii, sep="_")
     fn <- paste(outdir, "1_comb_", comb, "_motif.txt", sep="")
     x <- read.table(fn)$V1
     x
  })
  motifs <- unique(unlist(motifs))
    
  ###
  opfn <- paste(outdir, "2_treats_", ii, "_motif.txt", sep="") 
  write.table(motifs, opfn, row.names=F, quote=F, col.names=F)
  ## 
  cat(ii, length(motifs), "\n")
}

###################################################
### For each cell-type union of response motifs ###
###################################################

MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
for (i in 1:4){
   ##
   oneMCl <- MCls[i]
   res2 <- res%>%filter(MCls==oneMCl)
   ## 
   if( oneMCl!="Monocyte"){ 
      th0 <- quantile(abs(res2$beta),probs=0.9)
   }else{   
      th0 <- quantile(abs(res$beta),probs=0.9)
   }
   ## 
   motifs <- res2%>%
       filter(qval<0.1, abs(beta)>th0)%>%pull(gene)%>%unique()
   ##
   opfn <- paste(outdir, "3_MCls_", oneMCl, "_union.motif.txt", sep="")
   write.table(motifs, opfn, row.names=F, quote=F, col.names=F)

   ##
   cat(oneMCl, th0, length(motifs), "\n")
}


##################################
### define  response motif (2) ###
##################################

### consider direction, positive and negative

## rm(list=ls())


## outdir <- "./2_motif.activities.outs/"

## res <- read_rds("./2_motif.activities.outs/3_motif.diff.results.rds")

## ##
## ii <- 4
## oneMCl <- "Tcell"
## res2 <- res%>%filter(MCls==oneMCl)

## th0 <- quantile(abs(res2$beta),probs=0.9)

## ### CTRL
## res3 <- res2%>%filter(contrast=="LPS"|contrast=="PHA", qval<0.1, abs(beta)>th0, beta<0)
## ###
## opfn <- paste(outdir, "6_", ii, "_",  oneMCl, "_CTRL.motif.txt", sep="")
## write.table(unique(res3$gene), opfn, row.names=F, quote=F, col.names=F)
## length(unique(res3$gene))

## ### LPS-EtOH
## res3.1 <- res2%>%filter(contrast=="LPS", qval<0.1, abs(beta)>th0, beta>0)
## ##
## res3.2 <- res2%>%filter(contrast=="LPS-DEX", qval<0.1, abs(beta)>th0, beta<0)
## res3 <- rbind(res3.1, res3.2)
## opfn <- paste(outdir, "6_",  ii, "_", oneMCl, "_LPS-EtOH.motif.txt", sep="")
## write.table(unique(res3$gene), opfn, row.names=F, quote=F, col.names=F)
## length(unique(res3$gene))


## ### LPS-DEX
## res3 <- res2%>%filter(contrast=="LPS-DEX", qval<0.1, abs(beta)>th0, beta>0)
## ##
## opfn <- paste(outdir, "6_",  ii, "_", oneMCl, "_LPS-DEX.motif.txt", sep="")
## write.table(unique(res3$gene), opfn, row.names=F, quote=F, col.names=F)
## length(unique(res3$gene)) 


## ### PHA-EtOH
## res3.1 <- res2%>%filter(contrast=="PHA", qval<0.1, abs(beta)>th0, beta>0)
## res3.2 <- res2%>%filter(contrast=="PHA-DEX", qval<0.1, abs(beta)>th0, beta<0) 
## res3 <- rbind(res3.1, res3.2)
## ##
## opfn <- paste(outdir, "6_", ii, "_", oneMCl, "_PHA-EtOH.motif.txt", sep="")
## write.table(unique(res3$gene), opfn, row.names=F, quote=F, col.names=F)
## length(unique(res3$gene))


## ### PHA-DEX
## res3 <- res2%>%filter(contrast=="PHA-DEX", qval<0.1, abs(beta)>th0, beta>0)
## opfn <- paste(outdir, "6_", ii, "_", oneMCl, "_PHA-DEX.motif.txt", sep="")
## write.table(unique(res3$gene), opfn, row.names=F, quote=F, col.names=F)
## length(unique(res3$gene))



#################################
### define response motif (3) ###
#################################

### union of treatment motifs

res <- read_rds("./2_motif.activities.outs/3_motif.diff.results.rds")
##

outdir <- "./2_motif.activities.outs/"


MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
for (i in 1:4){
   ##
   oneMCl <- MCls[i]
   res2 <- res%>%filter(MCls==oneMCl)
   ## 
   if( oneMCl!="Monocyte"){ 
      th0 <- quantile(abs(res2$beta),probs=0.9)
   }else{   
      th0 <- quantile(abs(res$beta),probs=0.9)
   }
   ## 
   df2 <- res2%>%
       filter(qval<0.1, abs(beta)>th0)%>%
       dplyr::select(gene, motif)%>%distinct(.keep_all=T)
   ##
   opfn <- paste(outdir, "7.", i, "_", oneMCl, "_union.motif.txt", sep="")
   write.table(df2, opfn, row.names=F, quote=F, col.names=T)

   ##
   cat(oneMCl, nrow(res2),  nrow(df2), th0, length(unique(res2$gene)), "\n")
}



##############################################################
### identify motif with high activities in some celll type ###
##############################################################
## atac2 <- subset(atac2, subset=MCls!="DC")

## Y <- atac2@assays$chromvar@data
## meta <- atac2@meta.data%>%
##    mutate(treat=gsub(".*-ATAC-|_.*", "", NEW_BARCODE),
##           bti=paste(MCls, treat, SNG.BEST.GUESS, sep="_"))%>%
##    dplyr::select(NEW_BARCODE, bti) 


###
### Differential activities results
resMotif <- read_rds("./2_motif.activities.outs/3_motif.diff.results.rds")%>%
    mutate(conditions=paste(MCls, contrast, sep="_"))
### threshold
MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
ths <- sapply(1:4, function(i){
   ##
   oneMCl <- MCls[i]
   res <- resMotif
   res2 <- res%>%filter(MCls==oneMCl) 
   if( oneMCl!="Monocyte"){ 
      th0 <- quantile(abs(res2$beta),probs=0.9)
   }else{   
      th0 <- quantile(abs(res$beta),probs=0.9)
   }
   th0
})
names(ths) <- MCls

### motif id and name    
motif_DF <- res%>%dplyr::select(gene, motif)%>%distinct(gene, .keep_all=T)

###
### motif occupancy 
## atac <- read_rds("./2_motif.activities.outs/1_scATAC.motifActivities.rds")
X <- read_rds("./1.3_motif.outs/0_motif.rds")

###
###
anno <- read_rds("../2_Differential/2.2_compareRNAandATAC.outs/2_annot.ChIPseeker.rds")%>%
   as.data.frame()%>%
   mutate(peak=paste(gsub("chr", "", seqnames), start, end, sep="-"))
anno2 <- anno%>%dplyr::select(peak, geneId)

###
### DEG
fn <- "/nfs/rprdata/julong/SCAIP/analyses/SCAIP-B1-6_2020.03.23/6_DEG.CelltypeNew_output/Filter2/2_meta.rds"
resDG <- read_rds(fn)%>%mutate(conditions=paste(MCls, contrast, sep="_"))
resDG <- resDG%>%dplyr::filter(qval<0.1, abs(beta)>0.5)


###
### DVG
fn <- "/nfs/rprdata/julong/SCAIP/analyses/SCAIP-B1-6_2020.03.23/10_RNA.Variance_output/tmp9/3_phiNew.meta"
resDV <- read.table(fn, header=T)%>%mutate(conditions=paste(MCls, contrast, sep="_"))
resDV <- resDV%>%dplyr::filter(qval<0.1, abs(beta)>0.5)


###
#### DARs
fn <- "../2_Differential/1.3_DiffPeak.outs/3.0_DESeq_indi.results.rds"
resDA <- read_rds(fn)%>%mutate(conditions=paste(MCls, contrast, sep="_"))%>%as.data.frame()
resDA <- resDA%>%dplyr::filter(p.adjusted<0.1, abs(estimate)>0.5)


###
### summary number of motifs
conditions <- sort(unique(resMotif$conditions))
###
motif_DF <- map_dfr(1:length(conditions), function(i){
   ###
   condition <- conditions[i]
   oneMCl <- gsub("_.*", "", condition)
   ##
   cat(condition, "\n") 
   res2 <- resMotif%>%filter(conditions==condition, abs(beta)>ths[oneMCl])
   resDG2 <- resDG%>%filter(conditions==condition)
   resDV2 <- resDV%>%filter(conditions==condition) 
   resDA2 <- resDA%>%filter(conditions==condition)
    
   ##
   motif_DF2 <- map_dfr(1:nrow(res2), function(k){
      ##
      id <- res2[k,"gene"]
      motif.name <-  res2[k, "motif"]
      ### 
      peak_motif <- rownames(X)[X[,id]==1]
      peak_DA <- resDA2$gene 
      gene_motif <- anno2%>%filter(peak%in%peak_motif, peak%in%peak_DA)%>%pull(geneId)%>%unique()
      ###
      tmp <- data.frame(motif.id=id, motif=motif.name,
                        n_DEG=sum(gene_motif%in%resDG2$gene), n_DVG=sum(gene_motif%in%resDV2$gene))
      tmp
   })    
   ##
   motif_DF2$conditions <- condition
   motif_DF2 
})


###
###
x1 <- motif_DF%>%group_by(conditions)%>%
    summarise(nmotif=n(),.groups="drop")%>%as.data.frame()

x2 <- motif_DF%>%filter(n_DEG>0, n_DVG==0)%>%group_by(conditions)%>%
    summarise(nmotif_DEGonly=n(),.groups="drop")%>%as.data.frame()

x3 <- motif_DF%>%filter(n_DEG==0, n_DVG>0)%>%group_by(conditions)%>%
    summarise(nmotif_DVGonly=n(),.groups="drop")%>%as.data.frame()
 
x4 <- motif_DF%>%filter(n_DEG>0, n_DVG>0)%>%group_by(conditions)%>%
    summarise(nmotif_DGboth=n(),.groups="drop")%>%as.data.frame()

##
DF_summ <- x1%>%left_join(x2, by="conditions")%>%
    left_join(x3, by="conditions")%>%
    left_join(x4, by="conditions")%>%
    mutate(nmotif_DEGonly=ifelse(is.na(nmotif_DEGonly), 0, nmotif_DEGonly),
           nmotif_DVGonly=ifelse(is.na(nmotif_DVGonly), 0, nmotif_DVGonly),
           nmotif_DGboth=ifelse(is.na(nmotif_DGboth), 0, nmotif_DGboth))
##
## opfn <- "./2_motif.activities.outs/8_motifs.xlsx"
## openxlsx::write.xlsx(DF_summ, opfn, overwrite=T)
opfn <- "./2_motif.activities.outs/8_motifs.csv"
write.csv(DF_summ, opfn, row.names=F)
    



