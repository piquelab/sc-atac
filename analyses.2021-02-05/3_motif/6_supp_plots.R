#####
####
library(tidyverse)
library(ggplot2)
library(cowplot)
library(patchwork)
library(grid)
library(gridExtra)
library(ggExtra)
library(RColorBrewer)
library(ggsci)
library(viridis)
library(ComplexHeatmap)
library(circlize)
library(ggtext)
library(glue)


rm(list=ls())

###
outdir <- "./6_pub.outs/2_supp_plots/"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)





### FigS3_1
### dot plots of DAR enriched in TF

fn <-  "./1.2_motif.outs/3_motif.enrich.direction.rds"
enrich <- read_rds(fn)
drt2 <- c("0"="Down", "1"="Up")
enrich <- enrich%>%mutate(direction2=drt2[as.character(direction)],
   cluster=paste(contrast, MCls, direction2, sep="."))


cl <- paste( rep(c("LPS", "PHA", "LPS+DEX", "PHA+DEX"),each=4),
   rep(c("Bcell", "Monocyte", "NKcell", "Tcell"), times=4), sep=".")
cl <- paste(rep(cl,times=2), rep(c("Up","Down"), each=16), sep=".")

cluster2 <- 1:32
names(cluster2) <- cl


## ii <- paste(rep(c("A","B","C","D"),each=4),rep(1:4,times=4), sep="")
## cl2 <- paste(rep(c("X","Y"),each=16), rep(ii,times=2), sep=".")
## cluster2 <- setNames(cl2, cl)
## lab2 <- setNames(gsub("-","+",cl),cl2)
col1 <- c("LPS"="#fb9a99", "LPS+DEX"="#e31a1c", "PHA"="#a6cee3", "PHA+DEX"="#1f78b4")
col2 <- c("Bcell"="#4daf4a", "Monocyte"="#984ea3", "NKcell"="#aa4b56", "Tcell"="#ffaa00")
MCl_lab <- c("Bcell"="B cell", "Monocyte"="Monocyte", "NKcell"="NK cell", "Tcell"="T cell")

enrich <- enrich%>%
  mutate(cluster=gsub("-", "+", cluster),
         contrast=gsub("-", "+", contrast),
         MCl2=MCl_lab[as.character(MCls)],
         col.contrast=col1[contrast],
         col.MCls=col2[MCls],
         ClusterValue=as.numeric(cluster2[cluster]),
         ClusterNew=glue("<i style='color:{col.contrast}'>{contrast}.<i style='color:{col.MCls}'>{MCl2}.<i style='color:black'>{direction2}"),
         ClusterNew=fct_reorder(ClusterNew, ClusterValue))


##
##
## res <- read_rds("./2_motif.activities.outs/3_motif.diff.results.rds")
## topmotif2 <- res%>%dplyr::filter(qval<0.1,abs(beta)>1.41)%>%dplyr::pull(motif)%>%unique()

###
topmotif <- enrich%>%
   ## dplyr::filter(fold.enrichment>1.41, qvalue.fisher<0.1)%>%
   ## dplyr::pull(motif)%>%unique()
   group_by(cluster)%>%
   top_n(n=6, wt=fold.enrichment)%>%ungroup()%>%dplyr::pull(motif)%>%unique()

###

enrich3 <- enrich%>%
   dplyr::filter(motif%in%topmotif, fold.enrichment>1.41, qvalue.fisher<0.1)
###
p <- ggplot(enrich3, aes(x=ClusterNew, y=motif.name))+
   geom_point(aes(size=fold.enrichment, colour=qvalue.fisher))+
   scale_colour_gradient(name="FDR",
      low="blue", high="red", na.value=NA, trans="reverse", n.breaks=5,
      guide=guide_colourbar(order=1))+
   scale_size_binned("fold enrichment",
      guide=guide_bins(show.limits=T, axis=T,
          axis.show=arrow(length=unit(1.5,"mm"), ends="both"), order=2),
      n.breaks=4)+
   ## scale_x_discrete(labels=lab2)+
    
   theme_bw()+
   theme(axis.title=element_blank(),
         axis.text.x=element_markdown(angle=45, hjust=1, size=12),
         axis.text.y=element_text(size=12),
         legend.title=element_text(size=12),
         legend.text=element_text(size=12))
###
figfn <- paste(outdir, "FigS3_1_dotplot.pdf", sep="")
pdf(figfn, width=8.2, height=8.8)
print(p)
dev.off()





### FigS3_2
### dot plots of DEG enriched in TF

fn <- "./4_motif.enrich.outs/1.3_motif.DEG.direction.Allpeaks.rds"
enrich <- read_rds(fn)%>%
   mutate(MCls=gsub("_.*", "", comb),
          contrast=gsub(".*[le]_|_[UD].*", "", comb),
          direction=gsub(".*_", "", comb),
          cluster=paste(contrast, MCls, direction, sep="."))%>%as.data.frame()


cl <- paste( rep(c("LPS", "PHA", "LPS+DEX", "PHA+DEX"),each=4),
   rep(c("Bcell", "Monocyte", "NKcell", "Tcell"), times=4), sep=".")
cl <- paste(rep(cl,times=2), rep(c("Up","Down"), each=16), sep=".")

cluster2 <- 1:32
names(cluster2) <- cl


## ii <- paste(rep(c("A","B","C","D"),each=4),rep(1:4,times=4), sep="")
## cl2 <- paste(rep(c("X","Y"),each=16), rep(ii,times=2), sep=".")
## cluster2 <- setNames(cl2, cl)
## lab2 <- setNames(gsub("-","+",cl),cl2)
col1 <- c("LPS"="#fb9a99", "LPS+DEX"="#e31a1c", "PHA"="#a6cee3", "PHA+DEX"="#1f78b4")
col2 <- c("Bcell"="#4daf4a", "Monocyte"="#984ea3", "NKcell"="#aa4b56", "Tcell"="#ffaa00")
MCl_lab <- c("Bcell"="B cell", "Monocyte"="Monocyte", "NKcell"="NK cell", "Tcell"="T cell")


enrich <- enrich%>%
  mutate(cluster=gsub("-", "+", cluster),
         contrast=gsub("-", "+", contrast),
         col.contrast=col1[contrast],
         col.MCls=col2[MCls],
         MCl2=MCl_lab[as.character(MCls)],         
         ClusterValue=as.numeric(cluster2[cluster]),
         ClusterNew=glue("<i style='color:{col.contrast}'>{contrast}.<i style='color:{col.MCls}'>{MCl2}.<i style='color:black'>{direction}"),
         ClusterNew=fct_reorder(ClusterNew, ClusterValue))


##
## cluster motif
fn <- "./2_motif.activities.outs/Figure1.4_row_cluster.txt"
geneCL <- read.table(fn, header=T)
geneCL <- geneCL%>%mutate(cluster=paste("cluster", cluster, sep=""))

df <- data.frame(cluster=c("cluster1", "cluster2", "cluster3", "cluster4"),
                 cluster2=c("4", "1", "2", "3"))
geneCL <- geneCL%>%left_join(df, by="cluster")

geneCL <- geneCL[,2:3]
names(geneCL) <- c("motif.name", "motif.value")


## top 5 motifs
## enrich2 <- enrich%>%
##    dplyr::filter(fold.enrichment>1, qvalue.hyper<0.1)%>%
##    group_by(ClusterNew)%>%
##    top_n(n=5, wt=fold.enrichment)%>%ungroup()

## topmotif <- enrich2%>%dplyr::pull(motif.name)
## topmotif <- unique(topmotif)

enrich3 <- enrich%>%
    dplyr::filter(motif.name%in%geneCL$motif.name, fold.enrichment>1, qvalue.hyper<0.1) 

enrich3 <- enrich3%>%left_join(geneCL)%>%
    mutate(motif.name2=fct_reorder(motif.name, as.numeric(motif.value), .desc=T)) 

###
p <- ggplot(enrich3, aes(x=ClusterNew, y=motif.name2))+
   geom_point(aes(size=fold.enrichment, colour=qvalue.hyper))+
   scale_colour_gradient(name="FDR",
      low="blue", high="red", na.value=NA, trans="reverse", n.breaks=5,
      guide=guide_colourbar(order=1))+
   scale_size_binned("fold enrichment",
      guide=guide_bins(show.limits=T, axis=T,
      axis.show=arrow(length=unit(1.5,"mm"), ends="both"), order=2),
      n.breaks=4)+
   ## scale_x_discrete(labels=lab2)+
   theme_bw()+
   theme(axis.title=element_blank(),
         axis.text.x=element_markdown(angle=45, hjust=1, size=12),
         axis.text.y=element_text(size=12),
         legend.title=element_text(size=12),
         legend.text=element_text(size=12))
         
###
figfn <- paste0(outdir, "FigS3_2_DEG_dotplot.pdf")
pdf(figfn, width=8.2, height=8.8)
print(p)
dev.off()



#############################
### FigS3_3, point  plots ###
#############################


### assign values to combination
MCl2 <- c("Bcell", "Monocyte", "NKcell", "Tcell")
contrast2 <- c("LPS", "PHA", "LPS+DEX", "PHA+DEX")

comb2 <- paste(rep(MCl2, times=4), rep(contrast2, each=4), sep="_")
comb_val <- 1:16
names(comb_val) <- comb2 


### correlation between LFC on DEG and changes on motif activities
fn <- "./6_pub.outs/1_main_plots/2_plots_DEGvsmotif.rds"
plotDF <- read_rds(fn)
plotDF <- plotDF%>%dplyr::rename("beta.x"="beta", "beta.y"="LFC.RNA")
 
df <- plotDF%>%group_by(MCls, contrast)%>%
    summarise(rr=cor(beta.x, beta.y, method="pearson"), .groups="drop")
df <- df%>%mutate(comb=gsub("-", "+", paste(MCls, contrast, sep="_")), category="DEG")


###
### correlation between LFC on DVG and changes on motif activities
fn <- "./6_pub.outs/1_main_plots/3_plots_DVGvsmotif.rds"
plotDF <- read_rds(fn)
plotDF <- plotDF%>%dplyr::rename("beta.x"="beta", "beta.y"="LFC.RNA")
df2 <- plotDF%>%group_by(MCls, contrast)%>%
    summarise(rr=cor(beta.x, beta.y, method="pearson"), .groups="drop")
df2 <- df2%>%mutate(comb=gsub("-", "+", paste(MCls, contrast, sep="_")), category="DVG")



## col1 <- c("Bcell"="#4daf4a", "Monocyte"="#984ea3", "NKcell"="#aa4b56", "Tcell"="#ffaa00")
## col2 <- c("LPS"="#fb9a99", "LPS+DEX"="#e31a1c", "PHA"="#a6cee3", "PHA+DEX"="#1f78b4")

###
## dfcomb <- rbind(df, df2)
## dfcomb <- dfcomb%>%mutate(comb_value=as.numeric(comb_val[comb]), comb2=fct_reorder(comb, comb_value))
## ##
## p <- ggplot(dfcomb, aes(x=comb, y=rr, color=category, shape=category))+
##     geom_line(aes(group=category))+
##     geom_point()+
##     scale_color_manual(values=c("DEG"="red", "DVG"="blue"),
##        labels=c("DEG"="Gene expression", "DVG"="Gene variability"))+
##     scale_shape_manual(values=c("DEG"=2, "DVG"=5),
##        labels=c("DEG"="Gene expression", "DVG"="Gene variability"))+
##     ylim(-1,1)+ylab("PCC")+
##     theme_bw()+
##     theme(legend.title=element_blank(),
##           axis.title.x=element_blank(),
##           axis.text.x=element_text(angle=45, hjust=1, size=9),
##           axis.title.y=element_text(size=10),
##           axis.text.y=element_text(size=10))

## ###
## figfn <- paste(outdir, "FigS3_3_point.png", sep="")
## png(figfn, width=550, height=380, res=100)
## print(p)
## dev.off()


p2 <- ggplot(dfcomb, aes(x=comb2, y=rr, color=category, shape=category))+
    ##geom_line(aes(group=category))+
    geom_point(shape=5, size=2.5)+
    geom_segment(aes(xend=comb, yend=0))+
    scale_color_manual(values=c("DEG"="red", "DVG"="blue"),
       labels=c("DEG"="Gene expression", "DVG"="Gene variability"))+
    ## scale_shape_manual(values=c("DEG"=2, "DVG"=5),
    ##    labels=c("DEG"="Gene expression", "DVG"="Gene variability"))+
    ylim(-1,1)+ylab("PCC")+
    theme_bw()+
    theme(legend.title=element_blank(),
          axis.title.x=element_blank(),
          axis.text.x=element_text(angle=45, hjust=1, size=9),
          axis.title.y=element_text(size=10),
          axis.text.y=element_text(size=12))

###
figfn <- paste(outdir, "FigS3_3_segment.pdf", sep="")
pdf(figfn, width=6.8, height=4.5)
print(p2)
dev.off()





############################################
### FigS3_3, heatmap showing correlation ###
############################################


###
### DEG 

fn <- "./6_pub.outs/1_main_plots/2_plots_DEGvsmotif.rds"
plotDF <- read_rds(fn)
plotDF <- plotDF%>%dplyr::rename("beta.x"="beta", "beta.y"="LFC.RNA")


comb <- sort(unique(plotDF$comb))
DFcorr <- plotDF%>% ##dplyr::filter(gene%in%topmotif)%>%
    group_by(comb)%>%
    summarise(rr=cor(beta.x, beta.y),
              pval=as.numeric(cor.test(beta.x, beta.y)$p.value), .groups="drop")%>%as.data.frame()

DFcorr2 <- DFcorr%>%
    mutate(MCls=gsub("_.*", "", comb), contrast=gsub(".*_", "", comb),
           is_sig=ifelse(pval<0.05, 1, NA),
           rr2=rr*is_sig)

mat <- DFcorr2%>%pivot_wider(id_cols=MCls, names_from=contrast, values_from=rr2)
## mat2 for plot with NA value
mat2 <- as.matrix(mat[,-1])
colnames(mat2) <- c("LPS", "LPS+DEX", "PHA", "PHA+DEX")
rownames(mat2) <- mat$MCls

mat2 <- mat2[, c("LPS", "PHA", "LPS+DEX", "PHA+DEX")]


### mat3 for adding text in plots
mat3 <- DFcorr2%>%pivot_wider(id_cols=MCls, names_from=contrast, values_from=rr)
MCls <- mat3$MCls
mat3 <- as.matrix(mat3[,-1])
colnames(mat3) <- c("LPS", "LPS+DEX", "PHA", "PHA+DEX")
rownames(mat3) <- MCls
mat3 <- mat3[, c("LPS", "PHA", "LPS+DEX", "PHA+DEX")]


###
###
###
## num <- sort(as.numeric(mat_olap))
## mybreak <- c(quantile(num[1:460], probs=seq(0, 1, length.out=5)),
##              quantile(num[461:484], probs=seq(0,1,length.out=5)))
 
mycol <- colorRamp2(seq(0, 1, length.out=12), colorRampPalette(c("white", "red"))(12))

p1 <- Heatmap(mat2, name="PCC", na_col="grey90", 
    col=mycol, cluster_rows=F, cluster_columns=F,
    row_names_gp=gpar(fontsize=12), column_names_gp=gpar(fontsize=12),
    heatmap_legend_param=list(at=c(0, 0.25, 0.5, 0.75, 1),
        grid_width=grid::unit(0.38, "cm"), legend_height=grid::unit(5, "cm"),
        title_gp=gpar(fontsize=10), labels_gp=gpar(fontsize=10)),
    cell_fun=function(j, i, x, y, width, height, fill){
       grid.text(round(mat3[i,j],digits=3), x, y, gp=gpar(fontsize=10))
    })

## p1 <- draw(p1)
## p1 <- grid.grabExpr(p1)



figfn <- paste(outdir, "FigS3_3_1_heatmap_DEG.png", sep="")
png(figfn, height=480, width=520, res=120)
print(p1)
dev.off()




###########
### DVG ###
###########

fn <- "./6_pub.outs/1_main_plots/3_plots_DVGvsmotif.rds"
plotDF <- read_rds(fn)
plotDF <- plotDF%>%dplyr::rename("beta.x"="beta", "beta.y"="LFC.RNA")


comb <- sort(unique(plotDF$comb))
DFcorr <- plotDF%>% ##dplyr::filter(gene%in%topmotif)%>%
    group_by(comb)%>%
    summarise(rr=cor(beta.x, beta.y),
              pval=as.numeric(cor.test(beta.x, beta.y)$p.value), .groups="drop")%>%as.data.frame()

DFcorr2 <- DFcorr%>%
    mutate(MCls=gsub("_.*", "", comb), contrast=gsub(".*_", "", comb),
           is_sig=ifelse(pval<0.05, 1, NA),
           rr2=rr*is_sig)

mat <- DFcorr2%>%pivot_wider(id_cols=MCls, names_from=contrast, values_from=rr2)
## mat2 for plot with NA value
mat2 <- as.matrix(mat[,-1])
colnames(mat2) <- c("LPS", "LPS+DEX", "PHA", "PHA+DEX")
rownames(mat2) <- mat$MCls

mat2 <- mat2[, c("LPS", "PHA", "LPS+DEX", "PHA+DEX")]


### mat3 for adding text in plots
mat3 <- DFcorr2%>%pivot_wider(id_cols=MCls, names_from=contrast, values_from=rr)
MCls <- mat3$MCls
mat3 <- as.matrix(mat3[,-1])
colnames(mat3) <- c("LPS", "LPS+DEX", "PHA", "PHA+DEX")
rownames(mat3) <- MCls
mat3 <- mat3[, c("LPS", "PHA", "LPS+DEX", "PHA+DEX")]


###
###
###
## num <- sort(as.numeric(mat_olap))
## mybreak <- c(quantile(num[1:460], probs=seq(0, 1, length.out=5)),
##              quantile(num[461:484], probs=seq(0,1,length.out=5)))

mycol <- colorRamp2(seq(-1, 1,length.out=20), colorRampPalette(rev(brewer.pal(n=7, name="RdBu")))(20))

p2 <- Heatmap(mat2, name="PCC", na_col="grey90", 
    col=mycol, cluster_rows=F, cluster_columns=F,
    row_names_gp=gpar(fontsize=12), column_names_gp=gpar(fontsize=12),
    heatmap_legend_param=list(at=seq(-1, 1, by=0.5),
        grid_width=grid::unit(0.38, "cm"), legend_height=grid::unit(5, "cm"),
        title_gp=gpar(fontsize=10), labels_gp=gpar(fontsize=10)),
    cell_fun=function(j, i, x, y, width, height, fill){
       grid.text(round(mat3[i,j],digits=3), x, y, gp=gpar(fontsize=10))
    })


figfn <- paste(outdir, "FigS3_3_2_heatmap_DVG.png", sep="")
png(figfn, height=480, width=520, res=120)
print(p2)
dev.off()












    


###
### FigS3_4-FigS3_6

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


fn <- "./6_pub.outs/1_main_plots/2_plots_DEGvsmotif.rds"
plotDF <- read_rds(fn)
plotDF <- plotDF%>%dplyr::rename("beta.x"="beta", "beta.y"="LFC.RNA")

#### add cluster for motif
fn <- "./2_motif.activities.outs/Figure1.4_row_cluster.txt"
geneCL <- read.table(fn, header=T)
geneCL <- geneCL%>%mutate(cluster=paste("cluster", cluster, sep=""))

df <- data.frame(cluster=c("cluster1", "cluster2", "cluster3", "cluster4"),
                 cluster2=c("cluster4", "cluster1", "cluster2", "cluster3"))
geneCL <- geneCL%>%left_join(df, by="cluster")
 
plotDF <- plotDF%>%left_join(geneCL[,2:3], by="motif")



### 
col_MCls <- c("Bcell"="#4daf4a", "Monocyte"="#984ea3", "NKcell"="#aa4b56", "Tcell"="#ffaa00","DC"="#828282")
MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
MCl_lab <- c("Bcell"="B cell", "Monocyte"="Monocyte", "NKcell"="NK cell", "Tcell"="T cell")

###
i0 <- 4
for (i in c(1, 3, 4)){

## plot data
oneMCl <- MCls[i]
oneMCl2 <- MCl_lab[oneMCl]
    
plotDF2 <- plotDF%>%dplyr::filter(MCls==oneMCl) ##, gene%in%topmotif)    

cat(i, oneMCl, "\n")
    
    
xmin <- min(plotDF2$beta.x, na.rm=T)
xmax <- max(plotDF2$beta.y, na.rm=T)
xpos <- xmin+0.45*(xmax-xmin)

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

## map_dbl(data,~xFun(.x, a=0.3)),
## map_dbl(data,~yFun(.x, a=0.98)))%>%   

oneMCl2 <- paste(oneMCl2, "(gene expression)")    
p <- ggplot(plotDF2)+
    geom_point(aes(x=beta.x, y=beta.y, color=factor(cluster2)), size=0.8)+
    geom_text(data=anno_df2, aes(x=xpos, y=ypos, label=eq), colour="blue", size=4, parse=T)+
    scale_color_manual(values=c("cluster1"="#1b9e77", "cluster2"="#d95f02",
                                "cluster3"="#e7298a", "cluster4"="#7570b3"),
       labels=c("cluster1"="GR", "cluster2"="CEBP", "cluster3"="AP-1/NFKB", "cluster4"="STAT/IRF"),
       guide=guide_legend(override.aes=list(size=2)))+
    facet_wrap(~contrast, nrow=2, ncol=2, scales="fixed",
        labeller=as_labeller(c("LPS"="LPS", "LPS-DEX"="LPS+DEX", "PHA"="PHA", "PHA-DEX"="PHA+DEX")))+
    scale_x_continuous("TF activities changes", expand=expansion(mult=0.12))+
    scale_y_continuous("LFC on TF regulated genes expression", expand=expansion(mult=0.12))+
    ggtitle(oneMCl2)+
    theme_bw()+
    theme(
          strip.text=element_text(size=14),
          axis.title=element_text(size=12),
          axis.text=element_text(size=12),
          plot.title=element_text(hjust=0.5, size=14),
          legend.title=element_blank(),
          legend.text=element_text(size=12))        

p2 <- p+geom_smooth(data=plotDF2, aes(x=beta.x, y=beta.y),
   method="lm", formula=y~x, color="blue", size=0.5, se=F)

## slides
figfn <- paste0(outdir, "FigS3_", i0, "_",  oneMCl, "_motif_DEG.pdf")
pdf(figfn, width=6.2, height=5)  
print(p2)
dev.off()
     
i0 <- i0+1
cat(oneMCl, "\n")    
} ###



### 
### FigS3_7, FigS3_8 and FigS3_9 
### LFC on gene variability vs LFC on motif activities 

fn <- "./6_pub.outs/1_main_plots/3_plots_DVGvsmotif.rds"
plotDF <- read_rds(fn)
plotDF <- plotDF%>%dplyr::rename("beta.x"="beta", "beta.y"="LFC.RNA")

#### add cluster for motif
fn <- "./2_motif.activities.outs/Figure1.4_row_cluster.txt"
geneCL <- read.table(fn, header=T)
geneCL <- geneCL%>%mutate(cluster=paste("cluster", cluster, sep=""))

df <- data.frame(cluster=c("cluster1", "cluster2", "cluster3", "cluster4"),
                 cluster2=c("cluster4", "cluster1", "cluster2", "cluster3"))
geneCL <- geneCL%>%left_join(df, by="cluster")
 
plotDF <- plotDF%>%left_join(geneCL[,2:3], by="motif")



### 
col_MCls <- c("Bcell"="#4daf4a", "Monocyte"="#984ea3", "NKcell"="#aa4b56", "Tcell"="#ffaa00","DC"="#828282")
MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
MCl_lab <- c("Bcell"="B cell", "Monocyte"="Monocyte", "NKcell"="NK cell", "Tcell"="T cell")


pos_df <- data.frame(MCls=MCls, xpos=c(-0.5, 0.2, -0.2, -0.5), ypos=c(1.2, 0.35, 1, 0.9))


### plots
i0 <- 7
for (i in c(1, 3, 4)){
    
## plot data
oneMCl <- MCls[i]
oneMCl2 <- MCl_lab[oneMCl]    
plotDF2 <- plotDF%>%dplyr::filter(MCls==oneMCl) ##, gene%in%topmotif)    

    
pos_df2 <- pos_df%>%dplyr::filter(MCls==oneMCl)
xpos <- as.numeric(pos_df2$xpos)
ypos <- as.numeric(pos_df2$ypos)

    
anno_df2 <- plotDF2%>%
    group_by(contrast)%>%
    nest()%>%
    mutate(corr=map(data, ~cor.test((.x)$beta.x, (.x)$beta.y, method="pearson")),
          eq=map(corr,feq),
          r2=map_dbl(corr,~(.x)$estimate),
          xpos=xpos, ypos=ypos)%>%
   dplyr::select(-data,-corr)

oneMCl2 <- paste(oneMCl2, "(gene variability)")    
p <- ggplot(plotDF2)+
    geom_point(aes(x=beta.x, y=beta.y, color=cluster2), size=0.8)+
    geom_text(data=anno_df2, aes(x=xpos, y=ypos, label=eq), colour="blue", size=4, parse=T)+
    scale_color_manual(values=c("cluster1"="#1b9e77", "cluster2"="#d95f02",
                                "cluster3"="#e7298a", "cluster4"="#7570b3"),
       labels=c("cluster1"="GR", "cluster2"="CEBP", "cluster3"="AP-1/NFKB", "cluster4"="STAT/IRF"),                
       ## labels=c("cluster1"="pattern 1", "cluster2"="pattern 2", "cluster3"="pattern 3", "cluster4"="pattern 4"),
       guide=guide_legend(override.aes=list(size=2)))+
    facet_wrap(~contrast, nrow=2, ncol=2, scales="fixed",
        labeller=as_labeller(c("LPS"="LPS", "LPS-DEX"="LPS+DEX", "PHA"="PHA", "PHA-DEX"="PHA+DEX")))+
    scale_x_continuous("TF activities changes", expand=expansion(mult=0.18))+
    scale_y_continuous("LFC on TF regulated genes variability", expand=expansion(mult=0.18))+
    ggtitle(oneMCl2)+
    theme_bw()+
    theme(legend.title=element_blank(),
          legend.text=element_text(size=12),
          strip.text=element_text(size=14),
          axis.title=element_text(size=12),
          axis.text=element_text(size=12),
          plot.title=element_text(hjust=0.5, size=14))
    
p2 <- p+geom_smooth(data=plotDF2, aes(x=beta.x, y=beta.y),
   method="lm", formula=y~x, color="blue", size=0.5, se=F)

## slides
figfn <- paste0(outdir, "FigS3_", i0, "_",  oneMCl, "_motif_DVG.pdf")
pdf(figfn, width=6.2, height=5)  
print(p2)
dev.off()
     
i0 <- i0+1
cat(i0, oneMCl, "\n")    
    
    
} ###










################
#### Heatmap ###
################

### correlation between LFC on RNA and LFC on motif activities

comb <- sort(unique(plotDF$comb))
DFcorr <- plotDF%>% ##dplyr::filter(gene%in%topmotif)%>%
    group_by(comb)%>%
    summarise(rr=cor(beta.x, beta.y),
              pval=as.numeric(cor.test(beta.x, beta.y)$p.value), .groups="drop")%>%as.data.frame()

DFcorr2 <- DFcorr%>%
    mutate(MCls=gsub("_.*", "", comb), contrast=gsub(".*_", "", comb),
           is_sig=ifelse(pval<0.05, 1, NA),
           rr2=rr*is_sig)

mat <- DFcorr2%>%pivot_wider(id_cols=MCls, names_from=contrast, values_from=rr2)
## mat2 for plot with NA value
mat2 <- as.matrix(mat[,-1])
colnames(mat2) <- c("LPS", "LPS+DEX", "PHA", "PHA+DEX")
rownames(mat2) <- c("Bcell", "Monocyte", "NKcell", "Tcell")

mat2 <- mat2[, c("LPS", "PHA", "LPS+DEX", "PHA+DEX")]

### mat3 for adding text in plots
mat3 <- DFcorr2%>%pivot_wider(id_cols=MCls, names_from=contrast, values_from=rr)
mat3 <- as.matrix(mat3[,-1])
mat3 <- mat3[, c("LPS", "PHA", "LPS-DEX", "PHA-DEX")]
colnames(mat3) <- c("LPS", "PHA", "LPS+DEX", "PHA+DEX")

###
###
###
## num <- sort(as.numeric(mat_olap))
## mybreak <- c(quantile(num[1:460], probs=seq(0, 1, length.out=5)),
##              quantile(num[461:484], probs=seq(0,1,length.out=5)))
 
mycol <- colorRamp2(seq(0, 1, length.out=12), colorRampPalette(c("white", "red"))(12))

p <- Heatmap(mat2, name="PCC", na_col="grey90", 
    col=mycol, cluster_rows=F, cluster_columns=F,
    row_names_gp=gpar(fontsize=12), column_names_gp=gpar(fontsize=12),
    heatmap_legend_param=list(at=c(0, 0.25, 0.5, 0.75, 1),
        grid_width=grid::unit(0.38, "cm"), legend_height=grid::unit(4, "cm"),
        title_gp=gpar(fontsize=10), labels_gp=gpar(fontsize=10)),
    cell_fun=function(j, i, x, y, width, height, fill){
       grid.text(round(mat3[i,j],digits=3), x, y, gp=gpar(fontsize=10))
    })

figfn <- paste(outdir, "Figure1.2_topmotif_heatmap_DEG.png", sep="")
png(figfn, height=480, width=520, res=120)
print(p)
dev.off()





################
#### Heatmap ###
################

### correlation between LFC on RNA and LFC on motif activities
##
plotDF <- plotDF%>%dplyr::rename("beta.x"="beta", "beta.y"="LFC.RNA")


comb <- sort(unique(plotDF$comb))
DFcorr <- plotDF%>% ##dplyr::filter(gene%in%topmotif)%>%
    group_by(comb)%>%
    summarise(rr=cor(beta.x, beta.y),
              pval=as.numeric(cor.test(beta.x, beta.y)$p.value), .groups="drop")%>%as.data.frame()

DFcorr2 <- DFcorr%>%
    mutate(MCls=gsub("_.*", "", comb), contrast=gsub(".*_", "", comb),
           is_sig=ifelse(pval<0.05, 1, NA),
           rr2=rr*is_sig)

mat <- DFcorr2%>%pivot_wider(id_cols=MCls, names_from=contrast, values_from=rr2)

mat2 <- as.matrix(mat[,-1])
colnames(mat2) <- c("LPS", "LPS+DEX", "PHA", "PHA+DEX")
rownames(mat2) <- c("Bcell", "Monocyte", "NKcell", "Tcell")
mat2 <- mat2[, c("LPS", "PHA", "LPS+DEX", "PHA+DEX")]


mat3 <- DFcorr2%>%pivot_wider(id_cols=MCls, names_from=contrast, values_from=rr)
mat3 <- as.matrix(mat3[,-1])
colnames(mat3) <- gsub("-", "+", colnames(mat3))
mat3 <- mat3[,c("LPS", "PHA", "LPS+DEX", "PHA+DEX")]

###
###
###
## num <- sort(as.numeric(mat_olap))
## mybreak <- c(quantile(num[1:460], probs=seq(0, 1, length.out=5)),
##              quantile(num[461:484], probs=seq(0,1,length.out=5)))
 
mycol <- colorRamp2(seq(-1, 1,length.out=20), colorRampPalette(rev(brewer.pal(n=7, name="RdBu")))(20))

## mycol <- colorRamp2(seq(0, 1, length.out=12), colorRampPalette(c("white", "red"))(12))
 
p <- Heatmap(mat2, name="PCC", na_col="grey90",
    col=mycol, cluster_rows=F, cluster_columns=F,
    row_names_gp=gpar(fontsize=12), column_names_gp=gpar(fontsize=12),
    heatmap_legend_param=list(at=c(-1, -0.5, 0, 0.5, 1),
        grid_width=grid::unit(0.38, "cm"),
        legend_height=grid::unit(4, "cm"),
        labels_gp=gpar(fontsize=10),
        title_gp=gpar(fontsize=10)),
    cell_fun=function(j, i, x, y, width, height, fill){
       grid.text(round(mat3[i,j],digits=3), x, y, gp=gpar(fontsize=10))
    })

figfn <- paste(outdir, "Figure2.2_topmotif_heatmap_DVG.png", sep="")
png(figfn, height=480, width=520, res=120)
print(p)
dev.off()
