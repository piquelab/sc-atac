###
###
library(tidyverse)
library(ggplot2)
library(cowplot)
library(grid)
library(gridExtra)
library(ggExtra)
library(RColorBrewer)
library(ggsci)
library(viridis)
library(ComplexHeatmap)
library(circlize)

rm(list=ls())

###
outdir <- "./6_pub.outs/1_main_plots/"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)




### Figure-3A
### show heatmap of significant motif

res <- read_rds("./2_motif.activities.outs/3_motif.diff.results.rds")

## b0 <- quantile(abs(res$beta), probs=0.90)

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


condition3 <- c("B cell_LPS", "B cell_PHA", "Monocyte_LPS", "Monocyte_PHA",
                "NK cell_LPS", "NK cell_PHA", "T cell_LPS", "T cell_PHA",
               "B cell_LPS+DEX", "B cell_PHA+DEX", "Monocyte_LPS+DEX", "Monocyte_PHA+DEX",
                "NK cell_LPS+DEX", "NK cell_PHA+DEX", "T cell_LPS+DEX", "T cell_PHA+DEX")

## condition2 <- c("Bcell_LPS", "Bcell_PHA", "Bcell_LPS+DEX", "Bcell_PHA+DEX",
##                  "Monocyte_LPS", "Monocyte_PHA","Monocyte_LPS+DEX", "Monocyte_PHA+DEX",
##                  "NKcell_LPS", "NKcell_PHA", "NKcell_LPS+DEX", "NKcell_PHA+DEX",
##                  "Tcell_LPS", "Tcell_PHA", "Tcell_LPS+DEX", "Tcell_PHA+DEX")

## condition2 <- c("Bcell_LPS", "Monocyte_LPS", "NKcell_LPS", "Tcell_LPS",  
##                 "Bcell_PHA", "Monocyte_PHA", "NKcell_PHA", "Tcell_PHA",  
##                 "Bcell_LPS+DEX", "Monocyte_LPS+DEX", "NKcell_LPS+DEX", "Tcell_LPS+DEX",
##                 "Bcell_PHA+DEX", "Monocyte_PHA+DEX", "NKcell_PHA+DEX", "Tcell_PHA+DEX")
 
mat2 <- mat2[, condition2]
colnames(mat2) <- condition3


b <- as.vector(mat2)
b0 <- quantile(b[b<0], probs=seq(0, 1, length.out=49))
b1 <- 0
b2 <- quantile(b[b>0], probs=seq(0, 1, length.out=49))
breaks <- c(b0, b1, b2)

## breaks <- quantile(b, probs=seq(0, 1, length.out=100),na.rm=T)
 
col_fun <-  colorRamp2(breaks,
   colorRampPalette(rev(brewer.pal(n=7, name="RdBu")))(99))

column_ha <- HeatmapAnnotation(
   celltype=gsub("_.*", "", condition3),
   contrast=gsub("-", "+", gsub(".*_", "", condition3)),
   col=list(
      celltype=c("B cell"="#4daf4a", "Monocyte"="#984ea3",
                 "NK cell"="#aa4b56", "T cell"="#ffaa00"),
      contrast=c("LPS"="#fb9a99", "LPS+DEX"="#e31a1c",
                  "PHA"="#a6cee3", "PHA+DEX"="#1f78b4")),
   annotation_legend_param=list(celltype=list(labels_gp=gpar(fontsize=10), title_gp=gpar(fontsize=10),
      grid_width=grid::unit(0.5, "cm"), grid_height=grid::unit(0.6, "cm") ),
                                contrast=list(labels_gp=gpar(fontsize=10), title_gp=gpar(fontsize=10),
      grid_width=grid::unit(0.5, "cm"), grid_height=grid::unit(0.6, "cm"))),
   annotation_name_gp=gpar(fontsize=12))


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
   Pattern=list(labels=c("GR", "CEBP", "AP-1/NFKB", "STAT/IRF"),
                labels_gp=gpar(fontsize=10), title_gp=gpar(fontsize=10), grid_width=grid::unit(0.5, "cm"),
                grid_height=grid::unit(0.6, "cm"))),
   show_annotation_name=F, simple_anno_size=grid::unit(0.5, "cm")) 

                        

fig <- Heatmap(mat2, col=col_fun,
   cluster_rows=T, cluster_columns=F,
   show_row_dend=T,  show_column_dend=F,
   top_annotation=column_ha,
   right_annotation=row_ha,
   heatmap_legend_param=list(title="Diff motif",
      title_gp=gpar(fontsize=10),
      labels_gp=gpar(fontsize=10),
      legend_height=grid::unit(3.8, "cm"),
      grid_width=grid::unit(0.5, "cm")),
   show_row_names=T, show_column_names=T,
   ##row_names_side="left",
   column_names_gp=gpar(fontsize=10), column_names_rot=-45,
   row_names_gp=gpar(fontsize=9),
   use_raster=F, raster_device="png")
  
figfn <- paste(outdir, "Figure3.1_heatmap.motif.activities.png", sep="")
png(figfn, height=900, width=900, res=120)
set.seed(0)
fig <- draw(fig)
dev.off()




###
### scatter plots, LFC between motif activities and gene expression
### Figure3B


### motif differential results
resMotif <- read_rds("./2_motif.activities.outs/3_motif.diff.results.rds")%>%
    mutate(comb=paste(MCls, contrast, sep="_"))%>%as.data.frame()

topmotif <- resMotif%>%dplyr::filter(qval<0.1, abs(beta)>1.41)%>%dplyr::pull(gene)%>%unique()

resMotif2 <- resMotif%>%dplyr::filter(gene%in%topmotif)


### meta
fn <- "/nfs/rprdata/julong/SCAIP/analyses/SCAIP-B1-6_2020.03.23/6_DEG.CelltypeNew_output/Filter2/2_meta.rds"
DE_meta <- read_rds(fn)%>%
    mutate(comb=paste(MCls, contrast, sep="_"))

fn <- "./4.2_motif.enrich.outs/Motif_gene_DAR/Motif_genes.txt"
TF_gene <- read.table(fn, header=T)


#################
### plot data ###
#################
 

comb <- sort(unique(resMotif2$comb)) 
plotDF <- map_dfr(comb, function(ii){
   ###
   oneMCl <- gsub("_.*", "", ii) 
   res2 <- resMotif2%>%dplyr::filter(comb==ii)
   motifList <- unique(res2$gene) 
   DG <- DE_meta%>%dplyr::filter(qval<0.1, abs(beta)>0.5, MCls==oneMCl)%>%dplyr::pull(gene)%>%unique()
    
   LFC.RNA <- map_dfr(motifList, function(mm){
       ###
       gene2 <- TF_gene%>%dplyr::filter(MCls==oneMCl, motif_ID==mm)%>%pull(gene) 
       LFC <- DE_meta%>%dplyr::filter(comb==ii, gene%in%gene2, gene%in%DG)%>%
           dplyr::pull(beta)
       ###
       if ( length(LFC)>0){
          ## 
          tmp <- data.frame(LFC.RNA=median(LFC), ngene=length(LFC)) 
       }else{
          tmp <- data.frame(LFC.RNA=0, ngene=0) 
       }
       tmp
   })
   ###
   res2 <- cbind(res2, LFC.RNA) 
   res2
})


opfn <- paste(outdir, "2_plots_DEGvsmotif.rds", sep="")
write_rds(plotDF, file=opfn)

##
fn <- paste(outdir, "2_plots_DEGvsmotif.rds", sep="")
plotDF <- read_rds(fn)
plotDF <- plotDF%>%dplyr::rename("beta.x"="beta", "beta.y"="LFC.RNA")

## plotDF <- plotDF[,1:10]

###
#### add cluster for motif
fn <- "./2_motif.activities.outs/Figure1.4_row_cluster.txt"
geneCL <- read.table(fn, header=T)
geneCL <- geneCL%>%mutate(cluster=paste("cluster", cluster, sep=""))

df <- data.frame(cluster=c("cluster1", "cluster2", "cluster3", "cluster4"),
                 cluster2=c("cluster4", "cluster1", "cluster2", "cluster3"))
geneCL <- geneCL%>%left_join(df, by="cluster")
 
plotDF <- plotDF%>%left_join(geneCL[,2:3], by="motif")



#######################################
### scatter plots facet by contrast ###
#######################################
 
###
feq <- function(x){
   r <- round(as.numeric(x$estimate),digits=3)
   p <- x$p.value
   if(p<0.001) symb <- "***"
   if(p>=0.001 & p<0.01) symb <- "**"
   if (p>=0.01 & p<0.05) symb <- "*"
   if(p>0.05) symb <- "NS"
  
   eq <- bquote(italic(R)==.(r)~","~.(symb))
   eq 
}

##
xFun <- function(dx,a=0.5){
   min1 <- min(dx$beta.x, na.rm=T)
   max1 <- max(dx$beta.x, na.rm=T)
   R <- max1-min1
   xpos <- min1+a*R
}
##
yFun <- function(dx,a=0.8){
   min1 <- min(dx$beta.y, na.rm=T)
   max1 <- max(dx$beta.y, na.rm=T)
   R <- max1-min1
   ypos <- min1+a*R
}


### 
MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")

 
## plot data
oneMCl <- "Monocyte"   
plotDF2 <- plotDF%>%dplyr::filter(MCls==oneMCl) ##, gene%in%topmotif)    


xmin <- min(plotDF2$beta.x, na.rm=T)
xmax <- max(plotDF2$beta.y, na.rm=T)
xpos <- xmin+0.55*(xmax-xmin)

ymin <- min(plotDF2$beta.y, na.rm=T)
ymax <- max(plotDF2$beta.y, na.rm=T)
ypos <- ymin+0.98*(ymax-ymin)    
    
anno_df2 <- plotDF2%>%
    group_by(contrast)%>%
    nest()%>%
    mutate(corr=map(data, ~cor.test((.x)$beta.x, (.x)$beta.y, method="pearson")),
          eq=map(corr,feq),
          r2=map_dbl(corr,~(.x)$estimate),
          xpos=xpos, ypos=ypos)%>%
   dplyr::select(-data,-corr)


    
p <- ggplot(plotDF2)+
    geom_point(aes(x=beta.x, y=beta.y, color=factor(cluster2)), size=0.8)+
    geom_text(data=anno_df2, aes(x=xpos, y=ypos, label=eq), colour="blue", size=3, parse=T)+
    scale_color_manual(values=c("cluster1"="#1b9e77", "cluster2"="#d95f02",
                                "cluster3"="#e7298a", "cluster4"="#7570b3"),
       labels=c("cluster1"="pattern 1", "cluster2"="pattern 2", "cluster3"="pattern 3", "cluster4"="pattern 4"),
       guide=guide_legend(override.aes=list(size=2)))+
    facet_wrap(~contrast, nrow=2, ncol=2, scales="fixed",
        labeller=as_labeller(c("LPS"="LPS", "LPS-DEX"="LPS+DEX", "PHA"="PHA", "PHA-DEX"="PHA+DEX")))+
    scale_x_continuous("TF activities changes", expand=expansion(mult=0.12))+
    scale_y_continuous("LFC on TF regulated genes expression", expand=expansion(mult=0.12))+
    ggtitle("Monocyte (gene expression)")+
    theme_bw()+
    theme(legend.position="none",
          ## legend.title=element_blank(),
          ## legend.text=element_text(size=12),
          strip.text=element_text(size=14),
          axis.title=element_text(size=12),
          axis.text=element_text(size=12),
          plot.title=element_text(hjust=0.5, size=14))
         ## panel.background=element_rect(color=col_MCls[oneMCl], size=1.5))
          ##axis.title=element_text(size=12),
          ##plot.title=element_text(hjust=0.5))

p2 <- p+geom_smooth(data=plotDF2, aes(x=beta.x, y=beta.y),
   method="lm", formula=y~x, color="blue", size=0.5, se=F)

## slides
figfn <- paste(outdir, "Figure3.2_",  oneMCl, "_motif_DEG.png", sep="")
png(filename=figfn, width=500, height=580, res=120)  
print(p2)
dev.off()




### Figure3C
### Differential gene variability results

fn <- "/nfs/rprdata/julong/SCAIP/analyses/SCAIP-B1-6_2020.03.23/10_RNA.Variance_output/tmp9/3_phiNew.meta"
DE_meta <- read.table(fn, header=T)%>%mutate(comb=paste(MCls, contrast, sep="_"))


#################
### plot data ###
#################

comb <- sort(unique(resMotif2$comb)) 
plotDF <- map_dfr(comb, function(ii){
   ###
   cat(ii, "\n") 
   oneMCl <- gsub("_.*", "", ii) 
   res2 <- resMotif2%>%dplyr::filter(comb==ii)
   motifList <- unique(res2$gene) 
   DG <- DE_meta%>%dplyr::filter(qval<0.1, abs(beta)>0.5, MCls==oneMCl)%>%dplyr::pull(gene)%>%unique()
    
   LFC.RNA <- map_dfr(motifList, function(mm){
       ###
       gene2 <- TF_gene%>%dplyr::filter(MCls==oneMCl, motif_ID==mm)%>%pull(gene) 
       LFC <- DE_meta%>%dplyr::filter(comb==ii, gene%in%gene2, gene%in%DG)%>%
           dplyr::pull(beta)
       ###
       if ( length(LFC)>0){
          ## 
          tmp <- data.frame(LFC.RNA=median(LFC), ngene=length(LFC)) 
       }else{
          tmp <- data.frame(LFC.RNA=0, ngene=0) 
       }
       tmp
   })
   ###
   res2 <- cbind(res2, LFC.RNA) 
   res2
})


opfn <- paste(outdir, "3_plots_DVGvsmotif.rds", sep="")
write_rds(plotDF, file=opfn)
  

fn <- paste(outdir, "3_plots_DVGvsmotif.rds", sep="")
plotDF <- read_rds(fn)
##
plotDF <- plotDF%>%dplyr::rename("beta.x"="beta", "beta.y"="LFC.RNA")


### add cluster label
fn <- "./2_motif.activities.outs/Figure1.4_row_cluster.txt"
geneCL <- read.table(fn, header=T)
geneCL <- geneCL%>%mutate(cluster=paste("cluster", cluster, sep=""))

df <- data.frame(cluster=c("cluster1", "cluster2", "cluster3", "cluster4"),
                 cluster2=c("cluster4", "cluster1", "cluster2", "cluster3"))
geneCL <- geneCL%>%left_join(df, by="cluster")
 
plotDF <- plotDF%>%left_join(geneCL[,2:3], by="motif")



#######################################
### scatter plots facet by contrast ###
#######################################
 
MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
pos_df <- data.frame(MCls=MCls, xpos=c(-0.5, 0.2, -0.2, -0.5), ypos=c(1.2, 0.35, 1, 0.9))
 
## plot data
oneMCl <- "Monocyte"   
plotDF2 <- plotDF%>%dplyr::filter(MCls==oneMCl) ##, gene%in%topmotif)    

## xmin <- min(plotDF2$beta.x, na.rm=T)
## xmax <- max(plotDF2$beta.y, na.rm=T)
## xpos <- xmin+0.55*(xmax-xmin)

pos_df2 <- pos_df%>%dplyr::filter(MCls==oneMCl)
xpos <- as.numeric(pos_df2$xpos)
ypos <- as.numeric(pos_df2$ypos)
    
## ymin <- min(plotDF2$beta.y, na.rm=T)
## ymax <- max(plotDF2$beta.y, na.rm=T)
## ypos <- ymin+0.99*(ymax-ymin)    
    
anno_df2 <- plotDF2%>%
    group_by(contrast)%>%
    nest()%>%
    mutate(corr=map(data, ~cor.test((.x)$beta.x, (.x)$beta.y, method="pearson")),
          eq=map(corr,feq),
          r2=map_dbl(corr,~(.x)$estimate),
          xpos=xpos, ypos=ypos)%>%
   dplyr::select(-data,-corr)

    
p <- ggplot(plotDF2)+
    geom_point(aes(x=beta.x, y=beta.y, color=cluster2), size=0.8)+
    geom_text(data=anno_df2, aes(x=xpos, y=ypos, label=eq), colour="blue", size=3, parse=T)+
    scale_color_manual(values=c("cluster1"="#1b9e77", "cluster2"="#d95f02",
                                "cluster3"="#e7298a", "cluster4"="#7570b3"),
       labels=c("cluster1"="pattern 1", "cluster2"="pattern 2", "cluster3"="pattern 3", "cluster4"="pattern 4"),
       guide=guide_legend(override.aes=list(size=2)))+
    facet_wrap(~contrast, nrow=2, ncol=2, scales="fixed",
        labeller=as_labeller(c("LPS"="LPS", "LPS-DEX"="LPS+DEX", "PHA"="PHA", "PHA-DEX"="PHA+DEX")))+
    scale_x_continuous("TF activities changes", expand=expansion(mult=0.18))+
    scale_y_continuous("LFC on TF regulated genes variability", expand=expansion(mult=0.18))+
    ggtitle("Monocyte (gene variability)")+
    theme_bw()+
    theme(legend.position="none", ## legend.title=element_blank(),
          ## legend.text=element_text(size=12),
          strip.text=element_text(size=14),
          axis.title=element_text(size=12),
          axis.text=element_text(size=12),
          plot.title=element_text(hjust=0.5, size=14))
    
p2 <- p+geom_smooth(data=plotDF2, aes(x=beta.x, y=beta.y),
   method="lm", formula=y~x, color="blue", size=0.5, se=F)

## slides
figfn <- paste(outdir, "Figure3.3_", oneMCl, "_motif_DVG.png", sep="")
png(filename=figfn, width=500, height=580, res=120)  
print(p2)
dev.off()






