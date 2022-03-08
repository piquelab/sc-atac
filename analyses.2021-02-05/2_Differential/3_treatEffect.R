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

    

######################
#### scatter plots ###
######################

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

res <- read_rds("./1.2_DiffPeak.outs/2.0_DESeq.results.rds")%>%
   as.data.frame()%>%
   mutate(rn2=paste(MCls, gene,  sep="_"))%>%
   dplyr::rename("beta"="estimate")%>% 
   dplyr::filter(MCls!="DC")%>%drop_na(beta, p.adjusted) 

## DP <- res%>%drop_na(p.adjusted)%>%
##    dplyr::filter(p.adjusted<0.1, abs(estimate)>0.5)%>%
##    dplyr::pull(gene)
## DP <- as.character(unique(DP))



### (1), beta from LPS-EtOH vs CTRL against beta from LPS-DEX vs LPS-EtOH  
dfa <- res%>%filter(contrast=="LPS")    
dfb <- res%>%filter(contrast=="LPS-DEX")%>%dplyr::select(rn2, beta, p.value, p.adjusted)
       
df1 <- dfa%>%inner_join(dfb, by="rn2")

anno_df1 <- df1%>%group_by(MCls)%>%
   nest()%>%
   mutate(corr=map(data, ~cor.test((.x)$beta.x, (.x)$beta.y, method="pearson")),
          eq=map(corr,feq),
          r2=map_dbl(corr,~(.x)$estimate),
          xpos=map_dbl(data,~xFun(.x,a=0.7)),
          ypos=map_dbl(data,~yFun(.x,a=1)))%>%
   dplyr::select(-data,-corr)
     
fig1 <- ggplot(df1, aes(x=beta.x, y=beta.y))+
   rasterise(geom_point(size=0.3, color="grey50"),dpi=300)+ 
   geom_text(data=anno_df1, aes(x=xpos, y=ypos, label=eq), colour="blue", size=3, parse=T)+ 
   facet_wrap(~MCls, nrow=2, scales="free")+         
   scale_x_continuous("LPS effect on gene expression", expand=expansion(mult=0.1))+
   scale_y_continuous("LPS+DEX effect on gene expression", expand=expansion(mult=0.1))+
   theme_bw()+
   theme(strip.background=element_blank(),
         axis.title=element_text(size=12))
fig1 <- fig1+geom_smooth(method="lm",formula=y~x, size=0.5, se=F)
                           
figfn <- paste(outdir, "Figure2.1_LPS.png", sep="")
png(filename=figfn, width=500, height=500, pointsize=12, res=120)  
print(fig1)
dev.off()

### (2), beta from PHA-EtOH vs CTRL against beta from PHA-DEX vs PHA-EtOH
dfa <- res%>%filter(contrast=="PHA")    
dfb <- res%>%filter(contrast=="PHA-DEX")%>%dplyr::select(rn2, beta, p.value, p.adjusted)
       
df2 <- dfa%>%inner_join(dfb,by="rn2")

anno_df2 <- df2%>%group_by(MCls)%>%
    nest()%>%
    mutate(corr=map(data, ~cor.test((.x)$beta.x, (.x)$beta.y, method="pearson")),
          eq=map(corr,feq),
          r2=map_dbl(corr,~(.x)$estimate),
          xpos=map_dbl(data,~xFun(.x, a=0.7)),
          ypos=map_dbl(data,~yFun(.x, a=1)))%>%
   dplyr::select(-data,-corr)

fig2 <- ggplot(df2, aes(x=beta.x,y=beta.y))+
    rasterise(geom_point(size=0.3, color="grey50"),dpi=300)+ 
    geom_text(data=anno_df2, aes(x=xpos, y=ypos, label=eq), colour="blue", size=3, parse=T)+ 
    facet_wrap(~MCls, nrow=2, scales="free")+
    scale_x_continuous("PHA effect on gene expression", expand=expansion(mult=0.1))+
    scale_y_continuous("PHA+DEX effect on gene expression", expand=expansion(mult=0.1))+
    theme_bw()+
    theme(strip.background=element_blank(),
          axis.title=element_text(size=12))
fig2 <- fig2+geom_smooth(method="lm",formula=y~x, size=0.5, se=F)
                           
figfn <- paste(outdir, "Figure2.2_PHA.png", sep="")
png(filename=figfn, width=500, height=500, res=120)  
print(fig2)
dev.off()

###
###
figfn <- paste(outdir, "Figure2.3_comb.pdf", sep="")
pdf(figfn, width=10, height=5, pointsize=8)
print(plot_grid(fig1, fig2, nrow=1, ncol=2, labels="AUTO", label_fontface="plain"))
dev.off()
