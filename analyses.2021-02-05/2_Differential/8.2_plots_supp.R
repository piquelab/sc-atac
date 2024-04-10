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
library(qqman)
library(qvalue)
##
library(DESeq2)
library(biobroom)
library(ashr)
library(GenomicRanges)
library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(Signac)
library(clusterProfiler)
library(org.Hs.eg.db)
library(annotables)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(EnsDb.Hsapiens.v75)
library(ChIPseeker,lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
###
library(ggplot2)
library(cowplot) ##,lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
library(grid)
library(gridExtra)
library(ggExtra)
library(RColorBrewer)
library(ggsci)
library(RColorBrewer)
library(viridis)
library(ggrastr)
library(openxlsx)
##
library(ggtext, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
library(glue)



outdir <- "./8_pub.outs/2_supp_plots/"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)



###################
### 1. MA plots ###
###################

figfn <- paste(outdir, "FigS2_1_MA.png", sep="")
png(figfn, width=1000, height=1100, res=120)
par(mar=c(4,4,2,2),mgp=c(2,1,0))
x <- matrix(1:20, 5, 4, byrow=T)
layout(x)

fn <- "./1.3_DiffPeak.outs/3.0_DESeq_indi.results.rds"
res <- read_rds(fn)%>%
       mutate(color=ifelse((p.adjusted<0.1)&(!is.na(p.adjusted)), T, F))
 
MCls <- c("Bcell",  "Monocyte", "NKcell", "Tcell", "DC")
MCl2 <- c("B cell", "Monocyte", "NK cell", "T cell", "DC")
for (i in 1:length(MCls)){
    oneMCl <- MCls[i]
##1
   res2 <- res%>%dplyr::filter(MCls==oneMCl, contrast=="LPS")%>%
           dplyr::select(baseMean, estimate, color, p.value, p.adjusted)
   print(plotMA(res2[,1:3], colLine="NA", main="LPS", cex.main=1.5, font.main=1, cex.axis=1, cex.lab=1.2))
### LPS-DEX    
   res2 <- res%>%dplyr::filter(MCls==oneMCl, contrast=="LPS-DEX")%>%
          dplyr::select(baseMean, estimate, color, p.value, p.adjusted)
   print(plotMA(res2[,1:3], colLine="NA", main="LPS+DEX", cex.main=1.5, font.main=1, cex.axis=1, cex.lab=1.2))
### PHA
   res2 <- res%>%dplyr::filter(MCls==oneMCl, contrast=="PHA")%>%
           dplyr::select(baseMean, estimate, color, p.value, p.adjusted)
   print(plotMA(res2[,1:3], colLine="NA", main="PHA", cex.main=1.5, font.main=1, cex.axis=1, cex.lab=1.2))
### PHA-DEX
   res2 <- res%>%dplyr::filter(MCls==oneMCl, contrast=="PHA-DEX")%>%
           dplyr::select(baseMean, estimate, color, p.value, p.adjusted)
   print(plotMA(res2[,1:3], colLine="NA", main="PHA+DEX", cex.main=1.5, font.main=1, cex.axis=1, cex.lab=1.2))
    
   print(mtext(MCl2[i], side=4, line=0.9, cex=1.5, font=1)) 
}
dev.off()




###################
### 2, qq plots ###
###################

fn <- "./1.3_DiffPeak.outs/3.0_DESeq_indi.results.rds"
res <- read_rds(fn)%>%as.data.frame()%>%
   drop_na(p.value)%>%
   mutate(comb=paste(MCls, contrast, sep="_"))

comb <- sort(unique(res$comb))
dfNew <- map_dfr(comb, function(ii){
  res2 <- res%>%dplyr::filter(comb==ii)
  ngene <- nrow(res2)
  res2 <- res2%>%
     arrange(p.value)%>%
     mutate(observed=-log10(p.value), expected=-log10(ppoints(ngene)))
  res2
})


###
### facet by contrast
###
## dfNew$MCls <- gsub("DC", "z_DC", dfNew$MCls)
lab_treat <- c("LPS"="LPS", "LPS-DEX"="LPS+DEX", "PHA"="PHA", "PHA-DEX"="PHA+DEX")
lab_MCl <- c("Bcell"="B cell", "Monocyte"= "Monocyte", "NKcell"="NK cell",
          "Tcell"="T cell", "z_DC"="DC")
##
col_treat <- c("LPS"="#fb9a99", "LPS-DEX"="#e31a1c", "PHA"="#a6cee3", "PHA-DEX"="#1f78b4")
col_MCl <- c("Bcell"="#4daf4a", "Monocyte"="#984ea3", "NKcell"="#aa4b56", "Tcell"="#ffaa00", "z_DC"="#828282")

###
p1 <- ggplot(dfNew, aes(x=expected, y=observed, color=MCls))+
   ggrastr::rasterise(geom_point(size=0.3),dpi=300)+    
   geom_abline(colour="red")+
   scale_color_manual(values=col_MCl, labels=lab_MCl,
      guide=guide_legend(override.aes=list(size=2)))+
   facet_wrap(~contrast, scales="free",
      labeller=labeller(contrast=lab_treat), ncol=4)+
   xlab(bquote("Expected"~ -log[10]~"("~italic(plain(P))~")"))+
   ylab(bquote("Observed"~-log[10]~"("~italic(plain(P))~")"))+
   ###scale_y_continuous(expand=expansion(mult=c(0,0.3)))+
   theme_bw()+
   theme(legend.title=element_blank(),
         legend.background=element_blank(),
         legend.box.background=element_blank(),
         legend.key.size=grid::unit(0.8, "lines"),
         legend.text=element_text(size=10),
         axis.title=element_text(size=12),
         axis.text=element_text(size=12),
         strip.text.x=element_text(size=14))
###
###
figfn <- paste(outdir, "FigS2_2_treat_qq.pdf", sep="")
pdf(figfn, width=10, height=3.5)
print(p1)
dev.off()


###########################################################
### scatter plots of LFC between immune stimuli and DEX ###
###########################################################


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

MCl_label <- c("Bcell"="B cell", "Monocyte"="Monocyte", "NKcell"="NK cell", "Tcell"="T cell")


### (1), beta from LPS-EtOH vs CTRL against beta from LPS-DEX vs LPS-EtOH   
dfa <- res%>%dplyr::filter(contrast=="LPS")    
dfb <- res%>%dplyr::filter(contrast=="LPS-DEX")%>%dplyr::select(rn2, beta, p.value, p.adjusted)
       
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
   geom_text(data=anno_df1, aes(x=xpos, y=ypos, label=eq), colour="blue", size=4, parse=T)+ 
   facet_wrap(~MCls, nrow=2, scales="free", labeller=as_labeller(MCl_label))+         
   scale_x_continuous("LPS effect on DAR", expand=expansion(mult=0.1))+
   scale_y_continuous("LPS+DEX effect on DAR", expand=expansion(mult=0.1))+
   theme_bw()+
   theme(strip.background=element_blank(),
         strip.text=element_text(size=14),
         axis.title=element_text(size=12),
         axis.text=element_text(size=12))

fig1 <- fig1+geom_smooth(method="lm",formula=y~x, size=0.5, se=F)
                           


### (2), beta from PHA-EtOH vs CTRL against beta from PHA-DEX vs PHA-EtOH
dfa <- res%>%dplyr::filter(contrast=="PHA")    
dfb <- res%>%dplyr::filter(contrast=="PHA-DEX")%>%dplyr::select(rn2, beta, p.value, p.adjusted)
       
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
    geom_text(data=anno_df2, aes(x=xpos, y=ypos, label=eq), colour="blue", size=4, parse=T)+ 
    facet_wrap(~MCls, nrow=2, scales="free", labeller=as_labeller(MCl_label))+
    scale_x_continuous("PHA effect on DAR", expand=expansion(mult=0.1))+
    scale_y_continuous("PHA+DEX effect on DAR", expand=expansion(mult=0.1))+
    theme_bw()+
    theme(strip.background=element_blank(),
          strip.text=element_text(size=14),
          axis.title=element_text(size=12),
          axis.text=element_text(size=12))
fig2 <- fig2+geom_smooth(method="lm",formula=y~x, size=0.5, se=F)


###
###
figfn <- paste(outdir, "FigS2_3_scatter_DEX.pdf", sep="")
pdf(figfn, width=10, height=5, pointsize=8)
print(plot_grid(fig1, fig2, nrow=1, ncol=2, labels="AUTO", label_fontface="plain", label_size=14))
dev.off()




################################################
### enrichment analysis for nearby genes (2) ###
################################################

cg <- read_rds("./4_enrichment.outs/1_enrichGO.rds")

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
MCl_lab <- c("Bcell"="Bcell", "Monocyte"="Monocyte", "NKcell"="NK cell", "Tcell"="T cell")
 
cg <- cg%>%
  mutate(Cluster=gsub("-", "+", Cluster),
         contrast=gsub("-", "+", contrast),
         MCl2=MCl_lab[as.character(MCls)],
         col.contrast=col1[contrast],
         col.MCls=col2[MCls],
         ClusterValue=as.numeric(cluster2[Cluster]),
         ClusterNew=glue("<i style='color:{col.contrast}'>{contrast}.<i style='color:{col.MCls}'>{MCl2}.<i style='color:black'>{direction}"),
         ClusterNew=fct_reorder(ClusterNew, ClusterValue),
         maxGSSize=as.numeric(gsub("/.*", "", BgRatio)),
         ngene=as.numeric(gsub(".*/", "", GeneRatio)) )
 
cg2 <- cg%>%dplyr::filter(maxGSSize>10, maxGSSize<500, ngene>10, p.adjust<0.05)

fig1 <- enrichplot::dotplot(cg2, x="ClusterNew", showCategory=2)+
   scale_y_discrete(labels=function(y) str_wrap(y, width=50))+ 
   theme(axis.title=element_blank(),
         axis.text.x=element_markdown(angle=45, hjust=1,size=14),
         axis.text.y=element_text(size=14),
         legend.text=element_text(size=12),
         legend.title=element_text(size=12))

###
## figfn <- paste(outdir, "FigS2_4_enriched.png", sep="")
## png(figfn, width=1300, height=1200, res=120)
## print(fig1)
## dev.off()
 
###
figfn <- paste(outdir, "FigS2_4_GO_enriched.pdf", sep="")
pdf(figfn, width=13, height=12)
print(fig1)
dev.off()



################################
### KEGG enrichment analysis ###
################################

ck <- read_rds("./4_enrichment.outs/2_enrichKEGG.rds")


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
MCl_lab <- c("Bcell"="Bcell", "Monocyte"="Monocyte", "NKcell"="NK cell", "Tcell"="T cell")

ck <- ck%>%
  mutate(Cluster=gsub("-", "+", Cluster),
         contrast=gsub("-", "+", contrast),
         MCl2=MCl_lab[as.character(MCls)],
         col.contrast=col1[contrast],
         col.MCls=col2[MCls],
         ClusterValue=as.numeric(cluster2[Cluster]),
         ClusterNew=glue("<i style='color:{col.contrast}'>{contrast}.<i style='color:{col.MCls}'>{MCl2}.<i style='color:black'>{direction}"),
         ClusterNew=fct_reorder(ClusterNew, ClusterValue),
         maxGSSize=as.numeric(gsub("/.*", "", BgRatio)),
         ngene=as.numeric(gsub(".*/", "", GeneRatio)) )

ck2 <- ck%>%dplyr::filter(maxGSSize>10, maxGSSize<500, ngene>10, p.adjust<0.1)
 
fig2 <- enrichplot::dotplot(ck2, x="ClusterNew", showCategory=4)+
   scale_y_discrete(labels=function(y) str_wrap(y, width=50))+     
   theme(axis.title=element_blank(),
         axis.text.x=element_markdown(angle=45, hjust=1,size=14),
         axis.text.y=element_text(size=14),
         legend.title=element_text(size=12),
         legend.text=element_text(size=12))

##
figfn <- paste(outdir, "FigS2_4_KEGG_enriched.pdf", sep="")
pdf(figfn, width=13, height=12)
print(fig2)
dev.off()



#############################################
###  scatter plots of DEG or DVGs vs DARs ###
#############################################



fn <- "./1.3_DiffPeak.outs/3.0_DESeq_indi.results.rds"
resDP <- read_rds(fn)%>%drop_na(p.value)%>%
   mutate(comb=paste(MCls, contrast, sep="_"))%>%
   dplyr::rename("peak"="gene")%>%as.data.frame() 


fn <- "./2.2_compareRNAandATAC.outs/2_annot.ChIPseeker.rds"
peakAnno <- read_rds(fn)%>%as.data.frame()%>%
   mutate(peak=paste(gsub("chr", "", seqnames), start, end, sep="-")) 
peakAnno2 <- peakAnno%>%
    dplyr::select(peak, DTSS=distanceToTSS, geneId, SYMBOL)

                  ## flank_geneIds, flank_gene_distances, SYMBOL)



### test if DEG is DARs 
fn <- "/nfs/rprdata/julong/SCAIP/analyses/SCAIP-B1-6_2020.03.23/6_DEG.CelltypeNew_output/Filter2/2_meta.rds"
## fn <- "/nfs/rprdata/julong/SCAIP/analyses/SCAIP-B1-6_2020.03.23/6_DEG.CelltypeNew_output/Response_reviews/3.2_meta.ind.rds" 
resDE <- read_rds(fn)%>%drop_na(pval)%>%
   mutate(comb=paste(MCls, contrast, sep="_"), is_sig=ifelse(qval<0.1&abs(beta)>0.5, 1, 0))%>%
   as.data.frame() 
resDE2 <- resDE%>%dplyr::filter(qval<0.1, abs(beta)>0.5)
DEG <- unique(resDE2$gene)


### test if DVGs is DARs
fn <- "/nfs/rprdata/julong/SCAIP/analyses/SCAIP-B1-6_2020.03.23/10_RNA.Variance_output/tmp9/3_phiNew.meta"
resDV <- read.table(fn, header=T) %>%drop_na(pval)%>%
   mutate(comb=paste(MCls, contrast, sep="_"), is_sig=ifelse(qval<0.1&abs(beta)>0.5, 1, 0))
resDV2 <- resDV%>%dplyr::filter(qval<0.1, abs(beta)>0.5)
DVG <- unique(resDV2$gene)



###
### setting colors 

col1 <- c("LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
   "PHA"="#a6cee3", "PHA-DEX"="#1f78b4")
col2 <- c("Bcell"="#4daf4a", "Monocyte"="#984ea3",
   "NKcell"="#aa4b56", "Tcell"="#ffaa00")


## function
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

###
feq2 <- function(x){
   r <- round(as.numeric(x$estimate),digits=3)
   p <- x$p.value
   if(p<0.001) symb <- "***"
   if(p>=0.001 & p<0.01) symb <- "**"
   if (p>=0.01 & p<0.05) symb <- "*"
   if(p>0.05) symb <- "NS"
  
   eq <- bquote(italic(rho)==.(r)~","~.(symb))
   eq 
}



##
xFun <- function(dx,a=0.5){
   min1 <- 0
   max1 <- max(dx$npeak, na.rm=T)
   R <- max1-min1
   xpos <- min1+a*R
}
##
yFun <- function(dx,a=0.8){
   min1 <- 0
   max1 <- max(dx$x, na.rm=T)
   R <- max1-min1
   ypos <- min1+a*R
}



###
### plot data

resDP2 <- resDP%>%dplyr::filter(p.adjusted<0.1, abs(estimate)>0.5)
summDP <- resDP2%>%group_by(comb)%>%summarise(npeak=n(), .groups="drop")%>%ungroup()

summDE <- resDE2%>%group_by(comb)%>%summarise(x=n(), .groups="drop")%>%ungroup()

summDV <- resDV2%>%group_by(comb)%>%summarise(x=n(), .groups="drop")%>%ungroup()

###
plotDF <- summDP%>%inner_join(summDE, by="comb")%>%mutate(gr="DEG")
plotDF2 <- summDP%>%inner_join(summDV, by="comb")%>%mutate(gr="DVG")
plotDF <- rbind(plotDF, plotDF2)

##
plotDF <- plotDF%>%
   mutate(MCls=gsub("_.*", "", comb), contrast=gsub(".*_", "", comb))

anno_df2 <- plotDF%>%
    group_by(gr)%>%
    nest()%>%
    mutate(corr=map(data, ~cor.test((.x)$x, (.x)$npeak, method="spearman")),
          eq=map(corr,feq2),
          r2=map_dbl(corr,~(.x)$estimate),
          xpos=map_dbl(data,~xFun(.x, a=0.35)),
          ypos=map_dbl(data,~yFun(.x, a=0.9)))%>%
   dplyr::select(-data,-corr)


  
p <- ggplot(plotDF, aes(x=npeak, y=x))+
    geom_point(aes(color=contrast, shape=MCls))+
    geom_text(data=anno_df2, aes(x=xpos, y=ypos, label=eq), size=4, parse=T)+    
    scale_color_manual(values=col1,
       labels=c("LPS"="LPS", "LPS-DEX"="LPS+DEX", "PHA"="PHA", "PHA-DEX"="PHA+DEX"),
       guide=guide_legend(order=1, override.aes=list(shape=1)))+
    scale_shape_manual(values=c("Bcell"=0, "Monocyte"=1, "NKcell"=5, "Tcell"=2),
       labels=c("Bcell"="B cell", "Monocyte"="Monocyte", "NKcell"="NK cell", "Tcell"="T cell"),
       guide=guide_legend(order=2))+
    ylab("Number of DEGs (DVGs)")+
    xlab("Number of DARs")+
    facet_wrap(~gr, ncol=2, scales="free")+
    theme_bw()+
    theme(axis.title=element_text(size=12),
          axis.text=element_text(size=10),
          strip.text=element_text(size=14),
          legend.title=element_blank(),
          legend.background=element_blank(),
          legend.box.background=element_blank(),
          legend.key.size=grid::unit(0.8, "lines"),
          legend.text=element_text(size=10))


###
## figfn <- paste(outdir, "Figure0_DARs_scatter.png", sep="")
## png(figfn, width=580, height=300, res=120)
## p
## dev.off()

figfn <- paste(outdir, "FigS2_5_DARs_scatter.pdf", sep="")
pdf(figfn, width=6, height=2.8)
p
dev.off()




######################################################
### scatter plots LFC on ATAC vs LFC on DEG or DVG ###
######################################################

 
outdir <- "./8_pub.outs/2_supp/"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)


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
   min1 <- dx$min.x[1]
   max1 <- dx$max.x[1]
   R <- max1-min1
   xpos <- min1+a*R
}
##
yFun <- function(dx,a=0.8){
   min1 <- dx$min.y[1]
   max1 <- dx$max.y[1]
   R <- max1-min1
   ypos <- min1+a*R
}



fn <- "./1.3_DiffPeak.outs/3.0_DESeq_indi.results.rds"
resDP <- read_rds(fn)%>%drop_na(p.value)%>%
   mutate(comb=paste(MCls, contrast, sep="_"))%>%
   dplyr::rename("peak"="gene")%>%as.data.frame() 

## peakSel <- resDP%>%dplyr::filter(p.adjusted<0.1, abs(estimate)>0.5)%>%pull(peak)%>%unique()


###
### annotation 
fn <- "./2.2_compareRNAandATAC.outs/2_annot.ChIPseeker.rds"
peakAnno <- read_rds(fn)%>%as.data.frame()%>%
   mutate(peak=paste(gsub("chr", "", seqnames), start, end, sep="-")) 
peakAnno2 <- peakAnno%>%
    dplyr::select(peak, DTSS=distanceToTSS, geneId, SYMBOL)



###
### scatter plot of LFC of DEG vs DP

fn <- "/nfs/rprdata/julong/SCAIP/analyses/SCAIP-B1-6_2020.03.23/6_DEG.CelltypeNew_output/Filter2/2_meta.rds"
DEG <- read_rds(fn)%>%drop_na(pval)%>%
    dplyr::filter(abs(beta)>0.5, qval<0.1)%>%pull(gene)%>%unique()

resDG <- read_rds(fn)%>%drop_na(pval)%>%dplyr::filter(gene%in%DEG)%>%
    mutate(comb=paste(MCls, contrast, sep="_"), comb2=paste(comb, gene, sep="_"),
           is_sig=ifelse(qval<0.1&abs(beta)>0.5, 1, 0))


## fn <- "/nfs/rprdata/julong/SCAIP/analyses/SCAIP-B1-6_2020.03.23/6_DEG.CelltypeNew_output/Filter2/1_DESeq.results.rds"
## resDG <- read_rds(fn)%>%drop_na(p.value)%>%dplyr::filter(Batch2=="SCAIP6", gene%in%DEG)%>%
##     mutate(comb=paste(MCls, contrast, sep="_"), comb2=paste(comb, gene, sep="_"))%>%
##     as.data.frame()
## tmp <- resDE%>%mutate(comb2=paste(comb, gene, sep="_"))%>%dplyr::select(comb2, is_sig)
## resDG <- resDG%>%left_join(tmp, by="comb2")%>%dplyr::rename("beta"="estimate")


###
comb <- sort(unique(resDG$comb))    
plotdf <- map_dfr(comb, function(ii){
    ##
    cat(ii, "\n")
    x <- resDP%>%dplyr::filter(comb==ii)%>%
       left_join(peakAnno2, by="peak")%>%dplyr::filter(geneId%in%DEG) ##, peak%in%peakSel)

    ## geneDF <- lapply(1:nrow(x), function(k){
    ##     ##
    ##     genels <- unlist(str_split(x$flank_geneIds[k], ";"))
    ##     geneDTSS <- as.numeric(unlist(str_split(x$flank_gene_distances[k], ";")))
    ##     ##
    ##     DF <- data.frame(geneId=genels, DTSS=geneDTSS, peak=x$peak[k])
    ##     DF
    ##  })
    ##  geneDF <- do.call(rbind, geneDF)
        
    ##
    x2 <- x%>%
        dplyr::filter(abs(DTSS)<1e+05)%>%
       group_by(geneId)%>%summarise(beta.x=median(estimate, na.rm=T), .groups="drop")
       
    y <- resDG%>%dplyr::filter(comb==ii)%>%
       dplyr::select(gene, beta.y=beta, contrast, MCls, comb, is_sig)

    ##
    df <- x2%>%inner_join(y, by=c("geneId"="gene"))%>%as.data.frame()
    df
})

## opfn <- paste(outdir, "1_scatter_DEGvsDAR.txt", sep="")
## write.table(plotdf, file=opfn, row.names=F, col.names=T, quote=F, sep="\t")


df2 <- plotdf%>%group_by(MCls)%>%mutate(min.x=min(beta.x,na.rm=T), max.x=max(beta.x, na.rm=T))%>%ungroup()%>%
   group_by(contrast)%>%mutate(min.y=min(beta.y,na.rm=T), max.y=max(beta.y, na.rm=T))%>%ungroup() 
anno_df2 <- df2%>%
    group_by(contrast, MCls)%>%
    nest()%>%
    mutate(corr=map(data, ~cor.test((.x)$beta.x, (.x)$beta.y, method="pearson")),
          eq=map(corr,feq),
          r2=map_dbl(corr,~(.x)$estimate),
          xpos=0, ##map_dbl(data,~xFun(.x, a=0.7)),
          ypos=5.5)%>% ## map_dbl(data,~yFun(.x, a=1)))%>%
   dplyr::select(-data,-corr)

##
p <- ggplot(plotdf, aes(x=beta.x,y=beta.y))+
    rasterise(geom_point(aes(color=factor(is_sig)), size=0.3),dpi=300)+ 
    geom_text(data=anno_df2, aes(x=xpos, y=ypos, label=eq), colour="blue", size=3, parse=T)+
    scale_color_manual(values=c("1"="red", "0"="grey50"),
                       labels=c("1"="DEG", "0"="not DEG"),
                       guide=guide_legend(override.aes=list(size=3)))+
    facet_grid(contrast~MCls, scales="fixed",
       labeller=labeller(
          contrast=c("LPS"="LPS", "LPS-DEX"="LPS+DEX", "PHA"="PHA", "PHA-DEX"="PHA+DEX"),
          MCls=c("Bcell"="B cell", "Monocyte"="Monocyte", "NKcell"="NK cell", "Tcell"="T cell") ))+
    scale_x_continuous("LFC on chromatin accessibility", expand=expansion(mult=0.1))+
    scale_y_continuous("LFC on gene express changes", expand=expansion(mult=0.1))+
    theme_bw()+
    theme(axis.title=element_text(size=12),
          axis.text=element_text(size=12),
          strip.text=element_text(size=14),
          legend.title=element_blank(),
          legend.text=element_text(size=12))
##
p2 <- p+geom_smooth(method="lm",formula=y~x, size=0.5, se=F)
                           
###
###
figfn <- paste(outdir, "FigS2_6_DEGvsATAC.pdf", sep="")
pdf(figfn, width=7, height=6)  
print(p2)
dev.off()




###
### scatter plot, LFC of DVG vs DP

fn <- "/nfs/rprdata/julong/SCAIP/analyses/SCAIP-B1-6_2020.03.23/10_RNA.Variance_output/tmp9/3_phiNew.meta"
DVG <- read.table(fn, header=T)%>%drop_na(pval)%>%
    dplyr::filter(qval<0.1, abs(beta)>0.5)%>%pull(gene)%>%unique()
##
resDG <- read.table(fn, header=T) %>%drop_na(pval)%>%
   mutate(comb=paste(MCls, contrast, sep="_"),
          comb2=paste(comb, gene, sep="_"),
          is_sig=ifelse(qval<0.1&abs(beta)>0.5, 1, 0))%>%
   dplyr::filter(gene%in%DVG)

## fn <- "/nfs/rprdata/julong/SCAIP/analyses/SCAIP-B1-6_2020.03.23/10_RNA.Variance_output/tmp9/3_phiNew.results"
## resDG <- read.table(fn, header=T)%>%drop_na(pval)%>%dplyr::filter(batch=="SCAIP6", gene%in%DVG)%>%
##     mutate(comb=paste(MCls, contrast, sep="_"), comb2=paste(comb, gene, sep="_"))%>%
##     as.data.frame()
## tmp <- resDV%>%mutate(comb2=paste(comb, gene, sep="_"))%>%dplyr::select(comb2, is_sig)
## resDG <- resDG%>%left_join(tmp, by="comb2") ##%>%dplyr::rename("beta"="estimate")



comb <- sort(unique(resDG$comb))    
plotdf <- map_dfr(comb, function(ii){
    ##
    cat(ii, "\n")
    x <- resDP%>%dplyr::filter(comb==ii)%>%
       left_join(peakAnno2, by="peak")%>%dplyr::filter(geneId%in%DVG, peak%in%peakSel)
    ##
    x2 <- x%>%dplyr::filter(abs(DTSS)<1e+05)%>%
       group_by(geneId)%>%summarise(beta.x=median(estimate, na.rm=T), .groups="drop")
       
    y <- resDG%>%dplyr::filter(comb==ii)%>%
       dplyr::select(gene, beta.y=beta, contrast, MCls, comb, is_sig)
    
    ##
    df <- x2%>%inner_join(y, by=c("geneId"="gene"))%>%as.data.frame()
    df
})


df2 <- plotdf%>%group_by(MCls)%>%mutate(min.x=min(beta.x), max.x=max(beta.x))%>%ungroup()%>%
   group_by(contrast)%>%mutate(min.y=min(beta.y), max.y=max(beta.y))%>%ungroup() 
anno_df2 <- df2%>%
    group_by(contrast, MCls)%>%
    nest()%>%
    mutate(corr=map(data, ~cor.test((.x)$beta.x, (.x)$beta.y, method="spearman")),
          eq=map(corr,feq),
          r2=map_dbl(corr,~(.x)$estimate),
          xpos=0, ##map_dbl(data,~xFun(.x, a=0.6)),
          ypos=20)%>%  ##map_dbl(data,~yFun(.x, a=1)))%>%
   dplyr::select(-data,-corr)

x <- anno_df2%>%pivot_wider(id_cols=contrast, names_from=MCls, values_from=r2)%>%as.data.frame()

p <- ggplot(plotdf, aes(x=beta.x,y=beta.y))+
    rasterise(geom_point(aes(color=factor(is_sig)), size=0.3),dpi=300)+ 
    geom_text(data=anno_df2, aes(x=xpos, y=ypos, label=eq), colour="blue", size=3, parse=T)+
    scale_color_manual(values=c("1"="red", "0"="grey50"),
                       labels=c("1"="DVG", "0"="not DVG"),
                       guide=guide_legend(override.aes=list(size=3)))+
    facet_grid(contrast~MCls, scales="fixed",
       labeller=labeller(
          contrast=c("LPS"="LPS", "LPS-DEX"="LPS+DEX", "PHA"="PHA", "PHA-DEX"="PHA+DEX"),
            MCls=c("Bcell"="B cell", "Monocyte"="Monocyte", "NKcell"="NK cell", "Tcell"="T cell") ))+
    scale_x_continuous("LFC on chromatin accessibility", expand=expansion(mult=0.2))+  
    scale_y_continuous("LFC on gene variability", expand=expansion(mult=0.1))+
    theme_bw()+
    theme(axis.title=element_text(size=12),
          axis.text=element_text(size=12),
          strip.text=element_text(size=14),
          legend.title=element_blank(),
          legend.text=element_text(size=12))

p2 <- p+geom_smooth(method="lm",formula=y~x, size=0.5, se=F)
                           
## figfn <- paste(outdir, "Figure2.1_DVGvsATAC.png", sep="")
## png(filename=figfn, width=850, height=700, res=120)  
## print(p2)
## dev.off()

figfn <- paste(outdir, "FigS2_7_DVGvsATAC.pdf", sep="")
pdf(figfn, width=7, height=6)  
print(p2)
dev.off()


#################################
### Forest plots of DLDA_gene ###
#################################



MCls <- c("B cell", "Monocyte", "NK cell", "T cell")
contrast <- c("LPS", "LPS+DEX", "PHA", "PHA+DEX")
tmp2 <- data.frame(comb2=paste(rep(MCls, each=4), rep(contrast, times=4), sep="."))



col.treat <- c("LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
             "PHA"="#a6cee3", "PHA-DEX"="#1f78b4")
col.MCl <- c("Bcell"="#4daf4a", "Monocyte"="#984ea3",
             "NKcell"="#aa4b56", "Tcell"="#ffaa00")

MCl_label <- c("Bcell"="B cell", "Monocyte"="Monocyte", "NKcell"="NK cell", "Tcell"="T cell")

fn <- "./2.2_compareRNAandATAC.outs/3.3_enrich.pseudo.csv"
df <- read.csv(fn, header=T)
df3 <- data.frame(odds=df$odds,
   CI.low=df$lower, CI.high=df$upper,
   comb=gsub("-", "+", gsub("_", ".", df$comb)),
   MCls=df$cell, contrast=df$contrast, gene="DEG")

df3 <- df3%>%
   mutate(MCl2=MCl_label[as.character(MCls)],
          contrast2=gsub("-", "+", contrast),
          comb2=paste(MCl2, contrast2, sep="."))

df3 <- df3%>%full_join(tmp2, by="comb2")

p3 <- ggplot(df3, aes(x=odds, y=comb2))+
   geom_errorbarh(aes(xmax=CI.high, xmin=CI.low, colour=MCls),
       size=0.5, height=0.2)+
   geom_point(aes(colour=MCls), shape=19, size=1.5)+
   scale_colour_manual(values=col.MCl, labels=MCl_label)+
   geom_vline(aes(xintercept=1), size=0.25, linetype="dashed")+ 
   xlab("Odds ratio")+xlim(0, 30)+
   ## scale_y_discrete(labels=ylab2)+ 
   ggtitle("DLDA associated genes")+    
   theme_bw()+
   theme(plot.title=element_text(hjust=0.5, size=14),
         axis.title.y=element_blank(),
         axis.title.x=element_text(size=12),
         axis.text=element_text(size=12),
         legend.title=element_blank(),
         legend.text=element_text(size=12))


###
figfn <- paste(outdir, "FigS2_8_enrich_pseudo_forest.pdf", sep="")
pdf(figfn, width=6, height=4.5)
print(p3)
dev.off()




#############################
#### supplementary tables ###
#############################
library(data.table)
library(tidyverse)
outdir <- "./8_pub.outs/2_supp_plots/"

fn <- "./1.3_DiffPeak.outs/3.0_DESeq_indi.results.rds"
res <- read_rds(fn)%>%as.data.frame()%>%
   mutate(comb=paste(MCls, contrast, sep="_"))%>%
   dplyr::rename("peak"="gene") 

opfn <-paste(outdir, "TableS2_1_differential_summ.txt.gz", sep="")
fwrite(res, file=opfn, quote=F, sep=" ", na=NA)


###
fn <- paste(outdir, "TableS2_1_differential_summ.txt.gz", sep="")
x <- fread(fn, header=T, data.table=F, fill=T)



###
### binomial test

### Read data
fn <- "./1.3_DiffPeak.outs/3.0_DESeq_indi.results.rds"
res <- read_rds(fn)%>%as.data.frame()

##
res2 <- res%>%dplyr::filter(p.adjusted<0.1, abs(estimate)>0.5)%>%
    mutate(direction=sign(estimate), comb=paste(MCls, contrast, sep="_"))

combs <- sort(unique(res2$comb))
summ_df <- map_dfr(combs, function(ii){
    ##
    n1 <- res2%>%dplyr::filter(comb==ii, direction==-1)%>%pull(gene)%>%unique()%>%length()
    n2 <- res2%>%dplyr::filter(comb==ii, direction==1)%>%pull(gene)%>%unique()%>%length()
    ###
    ngene <- c(n1, n2)
    if(n1>n2){
       binom <- binom.test(ngene, 0.5, alternative="greater")
    }else{
       binom <- binom.test(ngene, 0.5, alternative="less") 
    }
    df2 <- data.frame(conditions=ii, neg=n1, npos=n2, nfold_neg_vs_pos=round(n1/n2, 2), 
                      nfold2=round(n2/n1,2), pval=format(binom$p.value, digits=4, scientific=T))
    df2
})

summ_df <- summ_df%>%mutate(conditions=gsub("-", "+", conditions))
###
opfn <- paste(outdir, "TableS2_3_binomial_test.xlsx", sep="")
write.xlsx(summ_df, file=opfn, overwrite=T)



