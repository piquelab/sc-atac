###
###
library(Matrix)
## library(MASS)
## library(scales)
library(tidyverse)
## library(parallel)
## library(data.table)
## library(future)
## library(purrr)
## library(furrr)
## library("BiocParallel")
## library(Rcpp)
## library(reshape)
## library(qqman)
## library(qvalue)

library(DESeq2)
## library(biobroom)
## library(ashr)
## library(GenomicRanges)
## library(Seurat)
## library(SeuratDisk)
## library(SeuratData)
## library(Signac)
## library(clusterProfiler)
## library(org.Hs.eg.db)
## library(annotables)
## library(TxDb.Hsapiens.UCSC.hg19.knownGene)
## library(EnsDb.Hsapiens.v75)
## library(ChIPseeker,lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
###
library(ggplot2)
library(cowplot) ##,lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
library(grid)
library(gridExtra)
library(ggExtra)
library(RColorBrewer)
library(ggsci)
library(RColorBrewer)
library(ComplexHeatmap)
library(viridis)
library(ggrastr)
library(openxlsx)
##
library(ggtext, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
library(glue)


#############################
### plots for publication ###
#############################

outdir <- "./8_pub.outs/1_main_plots/"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)


#######################################
### Figure 2 (A), barplots of #DARs ###
#######################################

fn <- "./1.3_DiffPeak.outs/3.0_DESeq_indi.results.rds"
res <- read_rds(fn)%>%as.data.frame()%>%mutate(direction=ifelse(estimate>0,1,0))
res2 <- res%>%dplyr::filter(p.adjusted<0.1, abs(estimate)>0.5)
sigs <- res2%>%group_by(MCls, contrast, direction)%>%
   summarise(ny=n(), .groups="drop")%>%
   mutate(ny2=ifelse(direction==0, -ny, ny))

breaks_value <- pretty(c(-12000, 15000), 8)

p <- ggplot(sigs, aes(x=MCls, y=ny2))+
   geom_bar(aes(fill=factor(MCls), alpha=factor(direction)), stat="identity")+
   scale_fill_manual(values=c("Bcell"="#4daf4a", "DC"="#828282",
      "Monocyte"="#984ea3", "NKcell"="#aa4b56", "Tcell"="#ffaa00"))+
   scale_alpha_manual(values=c("0"=0.5, "1"=1))+
   geom_hline(yintercept=0, color="grey60")+
   geom_text(aes(x=MCls, y=ny2, label=abs(ny2),
      vjust=ifelse(direction==1, -0.2, 1.2)), size=3)+
   scale_y_continuous("", breaks=breaks_value, limits=c(-13000,16000),
                      labels=abs(breaks_value))+
   scale_x_discrete(labels=c("Bcell"="B cell", "DC"="DC", "Monocyte"="Monocyte",
                             "NKcell"="NK cell", "Tcell"="T cell"))+ 
   facet_grid(~contrast,
      labeller=labeller(contrast=c("LPS"="LPS (CTRL)","LPS-DEX"="LPS+DEX (LPS)",
                                   "PHA"="PHA (CTRL)", "PHA-DEX"="PHA+DEX (PHA)")))+
   theme_bw()+
   theme(legend.position="none",
         strip.text.x=element_text(size=12),
         axis.title=element_blank(),
         axis.text.x=element_text(angle=-90, hjust=0, vjust=0.5, size=12),
         axis.text.y=element_text(size=12))

       
figfn <- "./8_pub.outs/1_main_plots/Figure2.1_barplot.png"
png(filename=figfn, width=850, height=400, res=120)
print(p)
dev.off()








############################################
### figure 2 (B), heatmap of correlation ###
############################################


res <- read_rds("./1.3_DiffPeak.outs/3.0_DESeq_indi.results.rds")%>%
   as.data.frame()%>%
   mutate(comb=paste(MCls, contrast, sep="_"))%>%
   dplyr::filter(MCls!="DC") 

DP <- res%>%drop_na(p.adjusted)%>%
   dplyr::filter(p.adjusted<0.1, abs(estimate)>0.5)%>%
   dplyr::pull(gene)%>%unique()%>%as.character()  ## 62323 peaks
## DP <- intersect(DP, peakSel) ## 47547 peaks

res2 <- res%>%dplyr::filter(gene%in%DP)

mat <- res2%>%pivot_wider(id_cols=gene, names_from=comb, values_from=estimate)
mat2 <- as.matrix(mat[,-1])

rnz <- rowSums(is.na(mat2))
mat2 <- mat2[rnz==0,]

colnames(mat2) <- gsub("-", "+", colnames(mat2))


###(1) heatmap
###mybreaks
 
## y <- as.vector(mat2)
## y0 <- y[abs(y)<1.73] #99% percent quantile(abs(y),probs=0.99)
## mybreaks <- c(min(y),quantile(y0,probs=seq(0,1,length.out=98)),max(y))
## names(mybreaks) <- NULL


col1 <- c("LPS"="#fb9a99", "LPS+DEX"="#e31a1c",
   "PHA"="#a6cee3", "PHA+DEX"="#1f78b4")
col2 <- c("B cell"="#4daf4a", "Monocyte"="#984ea3",
   "NK cell"="#aa4b56", "T cell"="#ffaa00")


### annotation
## x <- str_split(colnames(mat2), "_", simplify=T)
## anno_df <- data.frame("contrast"=x[,2],"celltype"=x[,1])

## col_ha <- HeatmapAnnotation(df=anno_df,
##     col=list(contrast=col1, celltype=col2),
##     annotation_name_gp=gpar(fontsize=12),
##     annotation_legend_param=list(
##   celltype=list(title_gp=gpar(fontsize=12), labels_gp=gpar(fontsize=10), title="celltype",
##                 grid_width=grid::unit(0.45, "cm"), grid_height=grid::unit(0.5, "cm")),        
##   contrast=list(title_gp=gpar(fontsize=12), labels_gp=gpar(fontsize=10), title="contrast",
##                 labels=c("LPS (ctrl)", "LPS+DEX (LPS)", "PHA (ctrl)", "PHA+DEX (PHA)"),
##                 grid_width=grid::unit(0.45, "cm"), grid_height=grid::unit(0.5, "cm"))))


## mycol <- colorRampPalette(rev(brewer.pal(n=7, name="RdBu")))(100)
## mycol2 <- colorRamp2(mybreaks, mycol)



###
### correlation of heatmap

Neworder <- c("Monocyte_LPS", "Monocyte_PHA", "Bcell_LPS", "Bcell_PHA",
              "NKcell_LPS", "NKcell_PHA", "Tcell_LPS", "Tcell_PHA",
              "Monocyte_LPS+DEX", "Monocyte_PHA+DEX", "Bcell_LPS+DEX", "Bcell_PHA+DEX",
              "NKcell_LPS+DEX", "NKcell_PHA+DEX", "Tcell_LPS+DEX", "Tcell_PHA+DEX")

Neworder2 <- c("Monocyte_LPS", "Monocyte_PHA", "B cell_LPS", "B cell_PHA",
              "NK cell_LPS", "NK cell_PHA", "T cell_LPS", "T cell_PHA",
              "Monocyte_LPS+DEX", "Monocyte_PHA+DEX", "B cell_LPS+DEX", "B cell_PHA+DEX",
              "NK cell_LPS+DEX", "NK cell_PHA+DEX", "T cell_LPS+DEX", "T cell_PHA+DEX")
 
corr <- cor(mat2)[Neworder, Neworder]
mycol <- colorRampPalette(rev(brewer.pal(n=7, name="RdBu")))(100) 
#mycol <- viridisLite::viridis(100)

colnames(corr) <- Neworder2
rownames(corr) <- Neworder2


x <- str_split(colnames(corr), "_", simplify=T)
df_col <- data.frame(celltype=x[,1], contrast=x[,2])
col_ha <- HeatmapAnnotation(df=df_col, col=list(contrast=col1, celltype=col2),
    annotation_name_gp=gpar(fontsize=12),
    annotation_legend_param=list(
  celltype=list(title_gp=gpar(fontsize=12), labels_gp=gpar(fontsize=10), title="celltype",
                grid_width=unit(0.45, "cm"), grid_height=unit(0.5, "cm")),
  contrast=list(title_gp=gpar(fontsize=12), labels_gp=gpar(fontsize=10), title="contrast",
                ## labels=c("LPS (CTRL)", "LPS+DEX (LPS)", "PHA (CTRL)", "PHA+DEX (PHA)"),
                grid_width=unit(0.45, "cm"), grid_height=unit(0.5, "cm"))),
  simple_anno_size=unit(0.4, "cm"))


row_ha <- rowAnnotation(df=df_col, col=list(contrast=col1, celltype=col2),
   show_annotation_name=F, show_legend=F,
   simple_anno_size=unit(0.4, "cm"))
 
p2 <- Heatmap(corr, col=mycol,
   cluster_rows=F, cluster_columns=F,
   show_row_names=F,
   show_column_names=T, column_names_gp=gpar(fontsize=12),
   column_names_rot=-90,   
   top_annotation=col_ha,
   left_annotation=row_ha,
   heatmap_legend_param=list(title="",
      title_gp=gpar(fontsize=12),
      labels_gp=gpar(fontsize=10),
      grid_width=grid::unit(0.35, "cm"),
      legend_height=grid::unit(6, "cm")),
   width=16*unit(5, "mm"), height=16*unit(5, "mm"),
   raster_device="png")

###
figfn <- "./8_pub.outs/1_main_plots/Figure2.2_heatmap.corr.png"
png(figfn, width=720, height=700,res=120)
set.seed(0)
p2 <- draw(p2)
dev.off()




##################################
### Figure 2(C), scatter plots ###
##################################

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
min1 <- min(dx$beta.x)
max2 <- max(dx$beta.x)
R <- max2-min1
xpos <- min1+a*R
}
##
yFun <- function(dx,a=0.8){
min1 <- min(dx$beta.y)
max2 <- max(dx$beta.y)
R <- max2-min1
ypos <- min1+a*R
}
  

    
### Read data
 
##load("./6_DEG.CelltypeNew_output/Filter2/Sigs.gene.DEG.RData")
outdir <- "./3_treatEffect.outs/"

res <- read_rds("./1.3_DiffPeak.outs/3.0_DESeq_indi.results.rds")%>%
   as.data.frame()%>%
   mutate(rn2=paste(MCls, gene,  sep="_"))%>%
   dplyr::rename("beta"="estimate")%>% 
   dplyr::filter(MCls!="DC")%>%drop_na(beta, p.adjusted) 

DP <- res%>%drop_na(p.adjusted)%>%
   dplyr::filter(p.adjusted<0.1, abs(beta)>0.5)%>%
   dplyr::pull(gene)
DP <- as.character(unique(DP))

res <- res%>%dplyr::filter(gene%in%DP)


### (2), beta from PHA-EtOH vs CTRL against beta from PHA-DEX vs PHA-EtOH
dfa <- res%>%filter(contrast=="PHA")    
dfb <- res%>%filter(contrast=="PHA-DEX")%>%dplyr::select(rn2, beta, p.value, p.adjusted)

###
###
df2 <- dfa%>%inner_join(dfb,by="rn2")
anno_df2 <- df2%>%group_by(MCls)%>%
    nest()%>%
    mutate(corr=map(data, ~cor.test((.x)$beta.x, (.x)$beta.y, method="pearson")),
          eq=map(corr,feq),
          r2=map_dbl(corr,~(.x)$estimate),
          xpos=map_dbl(data,~xFun(.x, a=0.7)),
          ypos=map_dbl(data,~yFun(.x, a=1)))%>%
   dplyr::select(-data,-corr)


#########################################
### T-cells PHA vs PHA-DEX as example ###
#########################################

###
plotDF <- df2%>%filter(MCls=="Tcell")
anno_plotDF <- anno_df2%>%filter(MCls=="Tcell")

p3 <- ggplot(plotDF, aes(x=beta.x,y=beta.y))+
    ## stat_density_2d(aes(fill=..level..), geom="polygon", contour=T)+ 
    rasterise(geom_point(size=0.3, color="grey50"),dpi=300)+ 
    geom_text(data=anno_plotDF, aes(x=xpos, y=ypos, label=eq), colour="blue", size=4, parse=T)+
    ## scale_fill_viridis_c(guide=guide_colorbar(barwidth=1, barheight=5))+
    scale_x_continuous("PHA effect on chromatin accessibility", expand=expansion(mult=0.1))+
    scale_y_continuous("PHA+DEX effect on chromatin accessibility", expand=expansion(mult=0.1))+
    geom_smooth(method="lm", formula=y~x, size=0.5, se=F)+
    ggtitle("T cell")+
    theme_minimal()+
    theme(##strip.background=element_blank(),
          legend.title=element_blank(),
          axis.title=element_text(size=12),
          axis.text=element_text(size=12),
          panel.background=element_rect(colour="#ffaa00",size=1.5),
          plot.title=element_text(hjust=0.5, size=14))

###                           
figfn <- "./8_pub.outs/1_main_plots/Figure2.3_Tcell.png"
png(filename=figfn, width=420, height=420, res=120)  
print(p3)
dev.off()



#############################
### Figure 2, (D) and (E) ###
#############################

### forest plots

col.treat <- c("LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
             "PHA"="#a6cee3", "PHA-DEX"="#1f78b4")
col.MCl <- c("Bcell"="#4daf4a", "Monocyte"="#984ea3",
             "NKcell"="#aa4b56", "Tcell"="#ffaa00")

###
### DEG if enriched in DARs

MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
contrast <- c("LPS", "LPS+DEX", "PHA", "PHA+DEX")
tmp2 <- data.frame(comb=paste(rep(MCls, each=4), rep(contrast, times=4), sep="."))



fn <- "./2.2_compareRNAandATAC.outs/3.1_enrich.DEG.csv"
df <- read.csv(fn, header=T)
df1 <- data.frame(odds=df$odds,
   CI.low=df$lower, CI.high=df$upper,
   comb=gsub("-", "+", gsub("_", ".", df$comb)),
   MCls=df$cell, contrast=df$contrast, gene="DEG")



df1 <- df1%>%full_join(tmp2, by="comb")


p1 <- ggplot(df1, aes(x=odds, y=comb))+
   geom_errorbarh(aes(xmax=CI.high, xmin=CI.low, colour=MCls),
       size=0.5, height=0.2)+
   geom_point(aes(colour=MCls), shape=19, size=1.5)+
   scale_colour_manual(values=col.MCl)+
   scale_y_discrete(
       labels=c("Bcell.LPS"="B cell.LPS", "Bcell.LPS+DEX"="B cell.LPS+DEX",
                "Bcell.PHA"="B cell.PHA", "Bcell.PHA+DEX"="B cell.PHA+DEX",
                "Monocyte.LPS"="Monocyte.LPS", "Monocyte.LPS+DEX"="Monocyte.LPS+DEX",
                "Monocyte.PHA"="Monocyte.PHA", "Monocyte.PHA+DEX"="Monocyte.PHA+DEX",
                "NKcell.LPS"="NK cell.LPS", "NKcell.LPS+DEX"="NK cell.LPS+DEX",
                "NKcell.PHA"="NK cell.PHA", "NKcell.PHA+DEX"="NK cell.PHA+DEX",
                "Tcell.LPS"="T cell.LPS", "Tcell.LPS+DEX"="T cell.LPS+DEX",
                "Tcell.PHA"="T cell.PHA", "Tcell.PHA+DEX"="T cell.PHA+DEX"))+ 
   geom_vline(aes(xintercept=1), size=0.25, linetype="dashed")+ 
   xlab("Odds ratio")+xlim(0, 10)+
   ## scale_y_discrete(labels=ylab2)+ 
   ggtitle("DEGs")+    
   theme_bw()+
   theme(plot.title=element_text(hjust=0.5, size=14),
         axis.title.y=element_blank(),
         axis.title.x=element_text(size=12),
         axis.text=element_text(size=10),
         legend.position="none")

###
fn <- "./2.2_compareRNAandATAC.outs/3.2_enrich.DVG.csv"
df <- read.csv(fn, header=T)
df2 <- data.frame(odds=df$odds,
   CI.low=df$lower, CI.high=df$upper,
   comb=gsub("-", "+", gsub("_", ".", df$comb)),
   MCls=df$cell, contrast=df$contrast)

p2 <- ggplot(df2, aes(x=odds, y=comb))+
   geom_errorbarh(aes(xmax=CI.high, xmin=CI.low, colour=MCls),
       size=0.5, height=0.2)+
   geom_point(aes(colour=MCls), shape=19, size=1.5)+
   scale_colour_manual(values=col.MCl)+
scale_y_discrete(
       labels=c("Bcell.LPS"="B cell.LPS", "Bcell.LPS+DEX"="B cell.LPS+DEX",
                "Bcell.PHA"="B cell.PHA", "Bcell.PHA+DEX"="B cell.PHA+DEX",
                "Monocyte.LPS"="Monocyte.LPS", "Monocyte.LPS+DEX"="Monocyte.LPS+DEX",
                "Monocyte.PHA"="Monocyte.PHA", "Monocyte.PHA+DEX"="Monocyte.PHA+DEX",
                "NKcell.LPS"="NK cell.LPS", "NKcell.LPS+DEX"="NK cell.LPS+DEX",
                "NKcell.PHA"="NK cell.PHA", "NKcell.PHA+DEX"="NK cell.PHA+DEX",
                "Tcell.LPS"="T cell.LPS", "Tcell.LPS+DEX"="T cell.LPS+DEX",
                "Tcell.PHA"="T cell.PHA", "Tcell.PHA+DEX"="T cell.PHA+DEX"))+     
   geom_vline(aes(xintercept=1), size=0.25, linetype="dashed")+
   xlab("Odds ratio")+xlim(0, 10)+
   ## scale_y_discrete(labels=ylab2)+   
   ggtitle("DVGs")+    
   theme_bw()+
   theme(plot.title=element_text(hjust=0.5, size=14),
         axis.title.y=element_blank(),
         axis.title.x=element_text(size=12),
         axis.text=element_text(size=10),
         legend.position="none")

figfn <- "./8_pub.outs/1_main_plots/Figure2.4_enrich_DEorDV_DAR.png"
png(figfn, width=850, height=420, res=120)
plot_grid(p1, p2, ncol=2)
dev.off()



##############
### poster ###
##############

rm(list=ls())

outdir2 <- "./8_pub.outs/1_main_plots/poster/"
if (!file.exists(outdir2)) dir.create(outdir2, showWarnings=F, recursive=T)

###
###
fn <- "./1.3_DiffPeak.outs/3.0_DESeq_indi.results.rds"
res <- read_rds(fn)%>%as.data.frame()%>%mutate(direction=ifelse(estimate>0,1,0))
res2 <- res%>%dplyr::filter(p.adjusted<0.1, abs(estimate)>0.5)
sigs <- res2%>%group_by(MCls, contrast, direction)%>%
   summarise(ny=n(), .groups="drop")%>%
   mutate(ny2=ifelse(direction==0, -ny, ny))

breaks_value <- pretty(c(-12000, 15000), 8)

p <- ggplot(sigs, aes(x=MCls, y=ny2))+
   geom_bar(aes(fill=factor(MCls), alpha=factor(direction)), stat="identity")+
   scale_fill_manual(values=c("Bcell"="#4daf4a", "DC"="#828282",
      "Monocyte"="#984ea3", "NKcell"="#aa4b56", "Tcell"="#ffaa00"))+
   scale_alpha_manual(values=c("0"=0.5, "1"=1))+
   geom_hline(yintercept=0, color="grey60")+
   geom_text(aes(x=MCls, y=ny2, label=abs(ny2),
      vjust=ifelse(direction==1, -0.2, 1.2)), size=3)+
   scale_y_continuous("", breaks=breaks_value, limits=c(-13000,16000),
                      labels=abs(breaks_value))+
   scale_x_discrete(labels=c("Bcell"="B cell", "DC"="DC", "Monocyte"="Monocyte",
                             "NKcell"="NK cell", "Tcell"="T cell"))+ 
   facet_grid(~contrast,
      labeller=labeller(contrast=c("LPS"="LPS","LPS-DEX"="LPS+DEX",
                                   "PHA"="PHA", "PHA-DEX"="PHA+DEX")))+
   theme_bw()+
   theme(legend.position="none",
         strip.text.x=element_text(size=14),
         axis.title=element_blank(),
         axis.text.x=element_text(angle=-45, hjust=0, vjust=0.5, size=10),
         axis.text.y=element_text(size=10),
         plot.margin=unit(c(5.5, 10, 5.5, 5.5), "points"))

       
figfn <- paste(outdir2, "Fig2.0_DAR.png", sep="")
png(filename=figfn, width=900, height=400, res=120)
print(p)
dev.off()








###
### Forest plot enrichmen analysis 

col.treat <- c("LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
             "PHA"="#a6cee3", "PHA-DEX"="#1f78b4")
col.MCl <- c("Bcell"="#4daf4a", "Monocyte"="#984ea3",
             "NKcell"="#aa4b56", "Tcell"="#ffaa00")


MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
contrast <- c("LPS", "LPS+DEX", "PHA", "PHA+DEX")
tmp2 <- data.frame(comb=paste(rep(MCls, each=4), rep(contrast, times=4), sep="."))

###
### DEG if enriched in DARs


fn <- "./2.2_compareRNAandATAC.outs/3.1_enrich.DEG.csv"
df <- read.csv(fn, header=T)
df1 <- data.frame(odds=log(df$odds),
   CI.low=log(df$lower), CI.high=log(df$upper),
   comb=gsub("-", "+", gsub("_", ".", df$comb)),
   MCls=df$cell, contrast=df$contrast, gene="DEG")

df1 <- df1%>%full_join(tmp2, by="comb")


p1 <- ggplot(df1, aes(x=odds, y=comb))+
   geom_errorbarh(aes(xmax=CI.high, xmin=CI.low, colour=MCls),
       linewidth=0.5, height=0.2)+
   geom_point(aes(colour=MCls), shape=19, size=1.5)+
   scale_colour_manual(values=col.MCl)+
   scale_y_discrete(
       labels=c("Bcell.LPS"="B cell.LPS", "Bcell.LPS+DEX"="B cell.LPS+DEX",
                "Bcell.PHA"="B cell.PHA", "Bcell.PHA+DEX"="B cell.PHA+DEX",
                "Monocyte.LPS"="Monocyte.LPS", "Monocyte.LPS+DEX"="Monocyte.LPS+DEX",
                "Monocyte.PHA"="Monocyte.PHA", "Monocyte.PHA+DEX"="Monocyte.PHA+DEX",
                "NKcell.LPS"="NK cell.LPS", "NKcell.LPS+DEX"="NK cell.LPS+DEX",
                "NKcell.PHA"="NK cell.PHA", "NKcell.PHA+DEX"="NK cell.PHA+DEX",
                "Tcell.LPS"="T cell.LPS", "Tcell.LPS+DEX"="T cell.LPS+DEX",
                "Tcell.PHA"="T cell.PHA", "Tcell.PHA+DEX"="T cell.PHA+DEX"))+ 
   geom_vline(aes(xintercept=0), size=0.25, linetype="dashed")+ 
   xlab("log odds ratio")+xlim(-1, 2.8)+
   ## scale_y_discrete(labels=ylab2)+ 
   ## ggtitle("DEGs")+    
   theme_bw()+
   theme(plot.title=element_text(hjust=0.5, size=14),
         axis.title.y=element_blank(),
         axis.title.x=element_text(size=12),
         axis.text=element_text(size=10),
         legend.position="none")

## figfn <- paste(outdir2, "Fig2.4_enrich_DE.png", sep="")
## png(figfn, width=420, height=420, res=120)
## p1
## dev.off()


###
fn <- "./2.2_compareRNAandATAC.outs/3.2_enrich.DVG.csv"
df <- read.csv(fn, header=T)
df2 <- data.frame(odds=log(df$odds),
   CI.low=log(df$lower), CI.high=log(df$upper),
   comb=gsub("-", "+", gsub("_", ".", df$comb)),
   MCls=df$cell, contrast=df$contrast)


df2 <- df2%>%full_join(tmp2, by="comb")

###
p2 <- ggplot(df2, aes(x=odds, y=comb))+
   geom_errorbarh(aes(xmax=CI.high, xmin=CI.low, colour=MCls),
       linewidth=0.5, height=0.2)+
   geom_point(aes(colour=MCls), shape=19, size=1.5)+
   scale_colour_manual(values=col.MCl)+
   scale_y_discrete(
       labels=c("Bcell.LPS"="B cell.LPS", "Bcell.LPS+DEX"="B cell.LPS+DEX",
                "Bcell.PHA"="B cell.PHA", "Bcell.PHA+DEX"="B cell.PHA+DEX",
                "Monocyte.LPS"="Monocyte.LPS", "Monocyte.LPS+DEX"="Monocyte.LPS+DEX",
                "Monocyte.PHA"="Monocyte.PHA", "Monocyte.PHA+DEX"="Monocyte.PHA+DEX",
                "NKcell.LPS"="NK cell.LPS", "NKcell.LPS+DEX"="NK cell.LPS+DEX",
                "NKcell.PHA"="NK cell.PHA", "NKcell.PHA+DEX"="NK cell.PHA+DEX",
                "Tcell.LPS"="T cell.LPS", "Tcell.LPS+DEX"="T cell.LPS+DEX",
                "Tcell.PHA"="T cell.PHA", "Tcell.PHA+DEX"="T cell.PHA+DEX"))+ 
   geom_vline(aes(xintercept=0), size=0.25, linetype="dashed")+ 
   xlab("log odds ratio")+xlim(-1, 2.8)+
   ## scale_y_discrete(labels=ylab2)+ 
   ## ggtitle("DVGs")+    
   theme_bw()+
   theme(plot.title=element_text(hjust=0.5, size=14),
         axis.title.y=element_blank(),
         axis.title.x=element_text(size=12),
         axis.text.x=element_text(size=10),
         axis.text.y=element_blank(),
         axis.ticks.y=element_blank(),
         legend.position="none")
 
figfn <- paste(outdir2, "Fig2.5_enrich_comb2.png", sep="")
png(figfn, width=700, height=380, res=120)
plot_grid(p1, p2, ncol=2, rel_widths=c(1, 0.65))
dev.off()



##################################
### bar plots of DEGs and DVGs ###
##################################

###
### DEG
fn <- "/nfs/rprdata/julong/SCAIP//analyses/SCAIP-B1-6_2020.03.23/6_DEG.CelltypeNew_output/Filter2/2_meta.rds"
res <- read_rds(fn)%>%filter(qval<0.1,abs(beta)>0.5)%>%drop_na(beta)

length(unique(res$gene))

sigs <- res%>%group_by(contrast, MCls)%>%
        summarise(ngene=n(),.groups="drop")%>%as.data.frame()
sigs <- sigs%>%mutate(contrast2=gsub("-", "+", contrast))


col2 <- c("Bcell"="#4daf4a", "Monocyte"="#984ea3", 
          "NKcell"="#aa4b56", "Tcell"="#ffaa00")

 
p0 <- ggplot(sigs, aes(x=contrast2, y=ngene))+
   geom_bar(aes(fill=MCls), stat="identity", position="dodge")+
   scale_fill_manual(values=col2)+
   theme_bw()+
   theme(legend.position="none",
         axis.title=element_blank(),
         axis.text=element_blank(),
         axis.ticks=element_blank())
 
figfn <- paste(outdir2, "Fig2.6_DEG.png", sep="")
ggsave(figfn, p0, width=400, height=180, units="px", dpi=120)
  

###
### DVG
fn <- "/nfs/rprdata/julong/SCAIP//analyses/SCAIP-B1-6_2020.03.23/10_RNA.Variance_output/tmp9/3_phiNew.meta"
res <- read.table(fn, header=T)%>%filter(qval<0.1,abs(beta)>0.5)%>%drop_na(beta)

length(unique(res$gene))

sigs <- res%>%group_by(contrast, MCls)%>%
        summarise(ngene=n(),.groups="drop")%>%as.data.frame()
sigs <- sigs%>%mutate(contrast2=gsub("-", "+", contrast))


col2 <- c("Bcell"="#4daf4a", "Monocyte"="#984ea3", 
          "NKcell"="#aa4b56", "Tcell"="#ffaa00")

 
p1 <- ggplot(sigs, aes(x=contrast2, y=ngene))+
   geom_bar(aes(fill=MCls), stat="identity", position="dodge")+
   scale_fill_manual(values=col2)+
   theme_bw()+
   theme(legend.position="none",
         axis.title=element_blank(),
         axis.text=element_blank(),
         axis.ticks=element_blank())

figfn <- paste(outdir2, "Fig2.7_DVG.png", sep="")
ggsave(figfn, p1, width=400, height=180, units="px", dpi=120)


