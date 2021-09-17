##
library(Matrix)
library(tidyverse)
library(DESeq2)
library(annotables)
##
library(ggplot2)
library(cowplot)
library(grid)
library(gridExtra)
library(ggExtra)
library(gtable)
library(ggsignif)
library(pheatmap)
library(ComplexHeatmap)
library(corrplot)
library(gtable)
library(RColorBrewer)
library(viridis)
library(ggrastr)

rm(list=ls())

##
outdir <- "./3_treatEffect.outs/"
if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=F)


#################################
### heatmap for 16 conditions ###
#################################

res <- read_rds("./1.2_DiffPeak.outs/2.0_DESeq.results.rds")%>%
   as.data.frame()%>%
   mutate(comb=paste(MCls, contrast, sep="_"))%>%
   dplyr::filter(MCls!="DC") 

DP <- res%>%drop_na(p.adjusted)%>%
   dplyr::filter(p.adjusted<0.1, abs(estimate)>0.5)%>%
   dplyr::pull(gene)
DP <- as.character(unique(DP))


###
comb <- unique(res$comb)
mat <- map_dfc(comb, function(ii){
##
  beta0 <- rep(NA, length(DP))
  names(beta0) <- DP
  d0 <- res%>%dplyr::filter(comb==ii, gene%in%DP)
  geneSel <- as.character(d0$gene)
  beta0[geneSel] <- d0$estimate
  beta0
})
mat <- as.matrix(mat)
rownames(mat) <- DP
colnames(mat) <- comb


col1 <- c("LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
   "PHA"="#a6cee3", "PHA-DEX"="#1f78b4")
col2 <- c("Bcell"="#4daf4a", "Monocyte"="#984ea3",
   "NKcell"="#aa4b56", "Tcell"="#ffaa00")


###
### heatmap-1
ii <- rowSums(is.na(mat))
mat2 <- mat[ii==0,]
y <- as.numeric(mat2)
y0 <- y[abs(y)<1.8] #99% percent quantile(abs(y),probs=0.99)
mybreaks <- c(min(y),quantile(y0,probs=seq(0,1,length.out=98)),max(y))
names(mybreaks) <- NULL



###
x <- str_split(comb, "_", simplify=T)
tmp_column <- data.frame(celltype=x[,1], treatment=x[,2])
rownames(tmp_column) <- comb
tmp_colors <- list(celltype=col2, treatment=col1) #brewer.pal(4,"Set1")
mycol <- colorRampPalette(rev(brewer.pal(n=7, name="RdBu")))(100)

#mycol <- viridisLite::viridis(100)
#mycol <- viridisLite::cividis(100, direction=1)
fig1 <- pheatmap(mat2, col=mycol, breaks=mybreaks, border_color="NA",
   cluster_rows=T, cluster_cols=T,
   annotation_col=tmp_column, annotation_colors=tmp_colors, annotation_legend=T,
   show_colnames=T, show_rownames=F, na_col="white",
   fontsize_row=12)

figfn <- "./3_treatEffect.outs/Figure1.1_heatmap.beta.png"
png(figfn, width=600, height=700,res=120)
print(fig1)
dev.off() 


###
### correlation heatmap
Neworder <- c("Monocyte_LPS", "Monocyte_PHA", "Bcell_LPS", "Bcell_PHA",
   "NKcell_LPS", "NKcell_PHA", "Tcell_LPS", "Tcell_PHA",
   "Monocyte_LPS-DEX", "Monocyte_PHA-DEX", "Bcell_LPS-DEX", "Bcell_PHA-DEX",
   "NKcell_LPS-DEX", "NKcell_PHA-DEX", "Tcell_LPS-DEX", "Tcell_PHA-DEX")
corr <- cor(mat2)[Neworder, Neworder]
mycol <- colorRampPalette(rev(brewer.pal(n=7, name="RdBu")))(100)
##

###
x <- str_split(colnames(corr), "_", simplify=T)
tmp_column <- data.frame(celltype=x[,1], treatment=x[,2])
rownames(tmp_column) <- colnames(corr)
tmp_colors <- list(celltype=col2, treatment=col1)

fig2 <- pheatmap(corr, col=mycol, scale="none", border_color="NA",
   cluster_rows=F, cluster_cols=F,
   annotation_col=tmp_column, annotation_colors=tmp_colors, annotation_legend =T,
   show_colnames=T, show_rownames=F, na_col="white")

###
figfn <- "./3_treatEffect.outs/Figure1.2_heatmap.corr.png"
png(figfn, width=600, height=600,res=120)
print(fig2)
dev.off() 

    
