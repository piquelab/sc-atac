###
library(tidyverse)
library(data.table)
library(qvalue)
##
library(Seurat)
library(Signac)
## library(SeuratDisk)
## library(SeuratData)
library(SeuratObject)

library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(viridis)
library(ComplexHeatmap)
library(circlize)
library(GGally)


rm(list=ls())

###
outdir <- "./3_summary.outs/"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)


###
### This script using for summary torus results from each motif results




################################
### summary number of motifs ###
################################

## MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")

## plotDF <- map_dfr(MCls, function(ii){
##     ##
##     fn <- paste("./torus_input/zzz_", ii, ".txt", sep="")
##     summ2 <- read.table(fn)
##     summ2
## })

fn <- "./torus_input/summary_comb/zzz_summ.txt"
plotDF <- read.table(fn)
names(plotDF) <- c("motif", "MCls", "num")

p <- ggplot(plotDF, aes(x=num))+
   geom_histogram(color="grey30", fill="grey", bins=50)+
   ## scale_color_manual(values=c("Bcell"="#4daf4a", "Monocyte"="#984ea3",
   ##     "NKcell"="#aa4b56", "Tcell"="#ffaa00"))+
   ##facet_wrap(~MCls, nrow=2, scales="free")+
   xlab("Number of SNPs for each motif")+
   ylab("Number of motifs")+
   theme_bw()+
   theme(legend.position="none",
         axis.title=element_text(size=12),
         axis.text=element_text(size=10))
         ##strip.text=element_text(size=12))
##   
figfn <- paste(outdir, "Figure0.1_numSNP_hist.png", sep="")
png(figfn, width=480, height=480, res=120)
p
dev.off()



##############################################################
### motif enrichment, 2 column                             ### 
### 1st column, 1 peaks and 0 otherwise                    ###
### 2rd column, 1 peaks also hit by motifs and 0 otherwise ###
##############################################################


motifList <- read.table("./Motif_file/motifList2020.txt")$V1

### motifs
DF <- map_dfr(motifList, function(ii){
###
  fn <- paste("./torus_output/comb/", ii, ".est", sep="")
  if ( file.exists(fn)&file.size(fn)>0){
     ### 
     est <- read.table(fn)%>%
        mutate(se=abs(V2-V3)/1.96, zval=V2/se, pval=pnorm(zval, lower.tail=F), motif=ii)%>%
        dplyr::rename(category=V1, estimate=V2, CI.lower=V3, CI.upper=V4)
     est2 <- est[3,] 
  }
})

###
## DF <- DF%>%mutate(FDR=qvalue(pval,pi0=1)$qvalue)%>%as.data.frame()

 
### motif name
fn <- "/nfs/rprdata/julong/sc-atac/analyses.2021-02-05/3_motif/2_motif.activities.outs/3_motif.diff.results.rds"
resDF <- read_rds(fn)
DF_motif <- resDF%>%distinct(gene, .keep_all=T)%>%
    dplyr::select(motif_ID=gene, motif_name=motif)

###
plotDF <- DF%>%left_join(DF_motif, by=c("motif"="motif_ID"))

###
### 18 motifs p<0.05 
plotDF2 <- plotDF%>%filter(pval<0.05) ##%>%slice_max(order_by=estimate, n=20)%>%as.data.frame()
plotDF2 <- plotDF2%>%mutate(motif_name2=fct_reorder(motif_name, estimate))
      
p <- ggplot(plotDF2, aes(x=estimate, y=motif_name2))+
   geom_errorbarh(aes(xmax=CI.upper, xmin=CI.lower), color="blue", size=0.5, height=0.2)+
   geom_point(shape=19, size=0.5, color="blue")+
   geom_vline(aes(xintercept=0), size=0.25, linetype="dashed")+
   scale_x_continuous("log odds ratio", breaks=seq(0, 10, 2), limits=c(-2,12))+
   theme_bw()+
   theme(legend.position="none",
         axis.title.x=element_text(size=10),
         axis.title.y=element_blank(),
         axis.text=element_text(size=10))

###
figfn <- paste(outdir, "Figure1_comb_top18_est.png", sep="")
png(figfn, width=380, height=500, res=120)
p
dev.off()


###
### extract peaks
## DF2 <- map_dfr(motifList, function(ii){
## ###
##   fn <- paste("./torus_output/comb/", ii, ".est", sep="")
##   if ( file.exists(fn)&file.size(fn)>0){
##      ### 
##      est <- read.table(fn)%>%
##         mutate(se=abs(V2-V3)/1.96, zval=V2/se, pval=pnorm(zval, lower.tail=F), motif=ii)%>%
##         dplyr::rename(category=V1, estimate=V2, CI_lower=V3, CI_upper=V4)
##      est2 <- est[2,] 
##   }
## })



###########################################
### summary annotation 2, single column ###
### 1-in peaks not by motifs            ###
### 2-in peaks hit by motfis            ###
### 0-otherwise                         ###
###########################################

motifList <- read.table("./Motif_file/motifList2020.txt")$V1

### motifs
DF <- map_dfr(motifList, function(ii){
###
  fn <- paste("./torus_output/comb2/", ii, ".est", sep="")
  if ( file.exists(fn)&file.size(fn)>0){
     ### 
     est <- read.table(fn)%>%
        mutate(se=abs(V2-V3)/1.96, zval=V2/se, pval=pnorm(zval, lower.tail=F), motif=ii)%>%
        dplyr::rename(category=V1, estimate=V2, CI.lower=V3, CI.upper=V4)
     est2 <- est[3,] 
  }
})


###
## DF <- DF%>%mutate(FDR=qvalue(pval,pi0=1)$qvalue)%>%as.data.frame()

  
### motif name
fn <- "/nfs/rprdata/julong/sc-atac/analyses.2021-02-05/3_motif/2_motif.activities.outs/3_motif.diff.results.rds"
resDF <- read_rds(fn)
DF_motif <- resDF%>%distinct(gene, .keep_all=T)%>%
    dplyr::select(motif_ID=gene, motif_name=motif)

###
plotDF <- DF%>%left_join(DF_motif, by=c("motif"="motif_ID"))


###
### 39 motifs p<0.05  
plotDF2 <- plotDF%>%filter(pval<0.05) ##%>%slice_max(order_by=estimate, n=20)%>%as.data.frame()
plotDF2 <- plotDF2%>%mutate(motif_name2=fct_reorder(motif_name, estimate))
      
p2 <- ggplot(plotDF2, aes(x=estimate, y=motif_name2))+
   geom_errorbarh(aes(xmax=CI.upper, xmin=CI.lower), color="blue", size=0.5, height=0.2)+
   geom_point(shape=19, size=0.5, color="blue")+
   geom_vline(aes(xintercept=0), size=0.25, linetype="dashed")+
   scale_x_continuous("log odds ratio", breaks=seq(0, 15, 3), limits=c(-2,15))+
   theme_bw()+
   theme(legend.position="none",
         ##axis.title.x=element_text(size=10),
         axis.title.y=element_blank(),
         ##axis.text.x=element_text(size=8),
         axis.text.y=element_text(size=8))

###
figfn <- paste(outdir, "Figure1.2_comb2_top39_est.png", sep="")
png(figfn, width=400, height=520, res=120)
p2
dev.off()



 

###################################
### Heatmap show motif activity ###
###################################

###
### 39 motifs pval<0.05, from above procedure
motif_top <- plotDF%>%filter(pval<0.05)%>%pull(motif_name)


### get motif activities
fn <- "/nfs/rprdata/julong/sc-atac/analyses.2021-02-05/3_motif/2_motif.activities.outs/1_scATAC.motifActivities.rds"
sc <- read_rds(fn)

meta <- sc@meta.data%>%filter(MCls=="Tcell")%>%
   mutate(treats=gsub("-", "+", gsub(".*-ATAC-|_.*", "", NEW_BARCODE)))%>%
   dplyr::select(NEW_BARCODE, MCls, treats)


Y <- sc@assays$chromvar@data
Y <- Y[, meta$NEW_BARCODE]
###
motif <- Motifs(sc)
motif_id <- rownames(Y)
motif_name <- ConvertMotifID(object=motif, id=motif_id)

rownames(Y) <- motif_name

## ### top motif list
## fn <- "./3_summary.outs/1_est.txt"
## DF_est <- read.table(fn, header=T)%>%drop_na(pval)%>%filter(category=="peaking.2")
 
## motif_top <- DF_est%>%filter(pval<0.05)%>%group_by(MCls)%>%slice_max(order_by=estimate, n=20)%>%as.data.frame()%>%
##     pull(motif_name)%>%unique()
##



Y2 <- Y[motif_top,]             
mat <- Y2 ##t(Y2)

anno_df <- meta%>%dplyr::select(celltype=MCls, treatment=treats)
ha <- HeatmapAnnotation(df=anno_df,
   col=list(celltype=c("Bcell"="#4daf4a", "Monocyte"="#984ea3", "NKcell"="#aa4b56", "Tcell"="#ffaa00"),
      treatment=c("CTRL"="#828282", "LPS"="#fb9a99", "LPS+DEX"="#e31a1c", "PHA"="#a6cee3", "PHA+DEX"="#1f78b4")),
   annotation_legend_param=list(
      celltype=list(grid_width=unit(0.3,"cm"),labels_gp=gpar(fontsize=8),
                    title_gp=gpar(fontsize=10), title="Cell type"),
      treatment=list(grid_width=unit(0.3,"cm"),labels_gp=gpar(fontsize=8),
                    title_gp=gpar(fontsize=10), title="Treats"))
   )


### setting color
b <- as.vector(mat)
bq <- quantile(b, probs=c(0.01, 0.99))
b1 <- quantile(b[b>=min(b)&b<bq[1]], probs=seq(0, 1, length.out=5))
b2 <- quantile(b[b>=bq[1]&b<bq[2]], probs=seq(0, 1, length.out=90))
b3 <- quantile(b[b>=bq[2]&b<max(b)], probs=seq(0, 1, length.out=5))
breaks <- c(b1, b2, b3)  
col_fun <- colorRamp2(breaks, colorRampPalette(brewer.pal(n=7, name="RdBu"))(100))


### main plots
fig <- Heatmap(mat, col=col_fun,
    cluster_rows=T, cluster_columns=F,
    show_row_dend=T, show_column_dend=F,
    show_row_names=T,
    show_column_names=F,
    row_names_gp=gpar(fontsize=8),
    top_annotation=ha,
    heatmap_legend_param=list(title="Motif activities",
       title_gp=gpar(fontsize=8),
       labels_gp=gpar(fontsize=8), grid_width=unit(0.3, "cm"), grid_height=unit(2, "cm")),
    use_raster=T, raster_device="png")
    

###
###
figfn <- paste(outdir, "Figure2_topmotif_heatmap.png", sep="")
png(figfn, height=500, width=1100, res=120)
set.seed(0)
fig <- draw(fig)
dev.off()



######################
### heatmap of LFC ###
######################

fn <- "/nfs/rprdata/julong/sc-atac/analyses.2021-02-05/3_motif/2_motif.activities.outs/3_motif.diff.results.rds"
res <- read_rds(fn)%>%
    mutate(comb=paste(MCls, contrast, sep="_"), is_sig=ifelse(qval<0.1, 1, 0))

matDF <- res%>%pivot_wider(id_cols=motif, names_from=comb, values_from=beta)

mat <- as.matrix(matDF[,-1])
rownames(mat) <- as.character(matDF$motif)

mat2 <- mat[as.character(motif_top),]

###
### is_sig
DF <- res%>%pivot_wider(id_cols=motif, names_from=comb, values_from=is_sig)
imat <- as.matrix(DF[,-1])
rownames(imat) <- as.character(DF$motif)
imat2 <- imat[as.character(motif_top),]


mat2 <- mat2*imat2


###
anno_df <- data.frame(rn=colnames(mat2))%>%
    mutate(celltype=gsub("_.*", "", rn), contrast=gsub("-", "+", gsub(".*_", "", rn)))%>%
    dplyr::select(-rn)
 
ha <- HeatmapAnnotation(df=anno_df,
   col=list(celltype=c("Bcell"="#4daf4a", "Monocyte"="#984ea3", "NKcell"="#aa4b56", "Tcell"="#ffaa00"),
      contrast=c( "LPS"="#fb9a99", "LPS+DEX"="#e31a1c", "PHA"="#a6cee3", "PHA+DEX"="#1f78b4")),
   annotation_legend_param=list(
      celltype=list(grid_width=unit(0.3,"cm"),labels_gp=gpar(fontsize=8),
                    title_gp=gpar(fontsize=10), title="cell-type"),
      contrast=list(grid_width=unit(0.3,"cm"),labels_gp=gpar(fontsize=8),
                    title_gp=gpar(fontsize=10), title="contrast"))
   )


### setting color
b <- as.vector(mat2)
b0 <- quantile(b[b<0], probs=seq(0, 1, length.out=49))
b1 <- 0
b2 <- quantile(b[b>0], probs=seq(0, 1, length.out=49))
breaks <- c(b0, b1, b2)

## bq <- quantile(b, probs=c(0.01, 0.99))
## b1 <- quantile(b[b>=min(b)&b<bq[1]], probs=seq(0, 1, length.out=5))
## b2 <- quantile(b[b>=bq[1]&b<bq[2]], probs=seq(0, 1, length.out=90))
## b3 <- quantile(b[b>=bq[2]&b<max(b)], probs=seq(0, 1, length.out=5))
## breaks <- c(b1, b2, b3)

col_fun <- colorRamp2(breaks, colorRampPalette(brewer.pal(n=7, name="RdBu"))(99))


### main plots
p2 <- Heatmap(mat2, col=col_fun,
    cluster_rows=T, cluster_columns=T,
    show_row_dend=T, show_column_dend=F,
    show_row_names=T,
    show_column_names=F,
    row_names_gp=gpar(fontsize=8),
    top_annotation=ha,
    heatmap_legend_param=list(title="LFC on motif",
       title_gp=gpar(fontsize=8),
       labels_gp=gpar(fontsize=8), grid_width=unit(0.3, "cm"), grid_height=unit(2, "cm")),
    use_raster=T, raster_device="png")
    

###
###
figfn <- paste(outdir, "Figure2.2_topmotif_LFC.png", sep="")
png(figfn, height=500, width=600, res=120)
set.seed(0)
p2 <- draw(p2)
dev.off()







##################################
### forest plots of odds ratio ###
##################################

## MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
## plotDF2 <- map_dfr(MCls, function(ii){
##    ##
##    fn <- paste("./torus_output/Peaks/Asthma_", ii, ".est", sep="")    
##    est <- read.table(fn)
##    est2 <- est[2,]%>%mutate(MCls=ii)%>%
##       dplyr::rename(category=V1, estimate=V2, CI.lower=V3, CI.upper=V4)
##    est2
## })    

## p <- ggplot(plotDF2, aes(x=estimate, y=MCls, color=factor(MCls)))+
##    geom_errorbarh(aes(xmax=CI.upper, xmin=CI.lower), size=0.5, height=0.2)+
##    geom_point(shape=19, size=0.5)+
##    geom_vline(aes(xintercept=0), size=0.25, linetype="dashed")+
##    scale_x_continuous("log odds ratio", breaks=seq(-2, 5, 2), limits=c(-2,5)) +
##    ylab("Cell-type active peaks")+ 
##    ## scale_y_discrete("SNP annotation", labels=motif_name2)+
##    scale_colour_manual(
##       values=c("Bcell"="#4daf4a", "Monocyte"="#984ea3",
##                "NKcell"="#aa4b56", "Tcell"="#ffaa00"))+
##    theme_bw()+
##    theme(legend.position="none",
##          axis.title=element_text(size=10),
##          axis.text=element_text(size=10))

## ###
## figfn <- paste(outdir, "Figure3_peaks_est.png", sep="")
## png(figfn, width=380, height=480, res=120)
## p
## dev.off()


## ###

## plotDF2 <- plotDF2%>%mutate(se=abs(CI.lower-estimate)/1.96, vi=(1/se)^2)
## b <- plotDF2$estimate
## vi <- plotDF2$vi
## bhat <- sum(vi*b)/sum(vi)
## se <- 1/sum(vi)



## ###############################################################
## ### summary combined peaks across cell-types for each motif ###
## ###############################################################

## ###
## ###
## motifList <- read.table("./Motif_file/motifList2020.txt")$V1

## ###
## DF_est <- map_dfr(motifList, function(ii){
## ###
##   fn <- paste("./torus_output/comb/", ii, ".est", sep="")
##   if ( file.exists(fn)&file.size(fn)>0){
##      ### 
##      est <- read.table(fn)%>%
##         mutate(se=abs(V2-V3)/1.96, zval=V2/se, pval=pnorm(abs(zval), lower.tail=F)*2, motif=ii)%>%
##         dplyr::rename(category=V1, estimate=V2, CI.lower=V3, CI.upper=V4)
##      est 
##   }
## })


## ###
## ### motif 
## fn <- "/nfs/rprdata/julong/sc-atac/analyses.2021-02-05/3_motif/2_motif.activities.outs/3_motif.diff.results.rds"
## resDF <- read_rds(fn)
## DF_motif <- resDF%>%distinct(gene, .keep_all=T)%>%
##     dplyr::select(motif_ID=gene, motif_name=motif)


## DF_est <- DF_est%>%left_join(DF_motif, by=c("motif"="motif_ID"))

## DF_est2 <- DF_est%>%filter(category=="peaking.2")%>%mutate(FDR=qvalue(pval)$qvalue)

## opfn <- "./3_summary.outs/2_comb_est.txt"
## write.table(DF_est2, opfn, quote=F, sep="\t", row.names=F, col.names=T)



## #########################################################
## ### summary motif with differential motif activities  ###
## #########################################################


## fn <- "./3_summary.outs/2_comb_est.txt"
## est <- read.table(fn, header=T)
## ##
## ##
## fn <- "/nfs/rprdata/julong/sc-atac/analyses.2021-02-05/3_motif/2_motif.activities.outs/3_motif.diff.results.rds"
## resDF <- read_rds(fn)


## ###
## MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
## mat <- NULL
## for (i in 1:4){
##    ##
##    oneMCl <- MCls[i]
##    res2 <- resDF%>%filter(MCls==oneMCl)
##    ## 
##    if( oneMCl!="Monocyte"){ 
##       th0 <- quantile(abs(res2$beta),probs=0.9)
##    }else{   
##       th0 <- quantile(abs(res2$beta),probs=0.9)
##    }
##    ##
##    ##
##    motifSel <- est%>%filter(pval<0.05)%>%pull(motif) 

##    df_pos <- res2%>%filter(qval<0.1, beta>th0, gene%in%motifSel)
##    tmp_pos <- df_pos%>%group_by(contrast)%>%summarise(nmotif=n(),.groups="drop")%>%ungroup()
##    tmp_pos <- tmp_pos%>%mutate(MCls=oneMCl, direction="pos") 
##    ### 
##    df_neg <- res2%>%filter(qval<0.1, beta<(-th0), gene%in%motifSel)
##    tmp_neg <- df_neg%>%group_by(contrast)%>%summarise(nmotif=n(),.groups="drop")
##    tmp_neg <- tmp_neg%>%mutate(MCls=oneMCl, direction="neg")
##    tmp <- rbind(tmp_pos, tmp_neg) 
##    ##
##    mat <- rbind(mat,tmp) 
## }
## ##


## rownames(mat) <- MCls
## colnames(mat) <- c("LPS", "LPS-DEX", "PHA", "PHA-DEX")







