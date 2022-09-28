#
library(Matrix)
library(MASS)
library(scales)
library(tidyverse)
library(parallel)
library(data.table)
##
library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(Signac)
#library(harmony)
library(annotables) 
library(org.Hs.eg.db)
###
library(ggplot2)
library(cowplot)
library(grid)
library(gridExtra)
library(ggExtra)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(viridis)
library(ggrastr)
theme_set(theme_grey())

rm(list=ls())

outdir <- "./5_dynamic.outs/"
if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=FALSE)


########################################
### motif activties dynamic analysis ###
###       6-29-2022, by JW           ###
########################################

sc <- read_rds("./2_motif.activities.outs/1_scATAC.motifActivities.rds")
motif <- Motifs(sc)
motif.names <- unlist(motif@motif.names)

MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")


###
### Differential motif
res <- read_rds("./2_motif.activities.outs/3_motif.diff.results.rds")
th0 <- quantile(abs(res$beta),probs=0.9)
## res2 <- res2%>%filter(qval<0.1, abs(beta)>th0)
## res2%>%group_by(MCls, contrast)%>%summarise(ny=n(),.groups="drop")


###
### extract data
for (oneMCl in MCls){

cat(oneMCl, "\n")
    
###    
sc2 <- subset(sc, subset=MCls==oneMCl)
X <- sc2@assays$chromvar@data
rn <- rownames(X)


meta <- sc2@meta.data%>%
       mutate(treats=gsub(".*-ATAC-|_.*", "", NEW_BARCODE))
treats <- unique(meta$treats)
ntreat <- length(treats)
ncell <- ncol(X)
np <- nrow(X) 

###
### gene variance across cells  
mu <- rowMeans(X)
Xc <- sweep(X, 1, mu, "-")
var <- rowSums(Xc*Xc)/(ncell-1)
names(var) <- rn


###
### calculate DLDA
contrast_ls <- list("LPS"=c("CTRL", "LPS"),
    "LPS-DEX"=c("LPS", "LPS-DEX"),
    "PHA"=c("CTRL", "PHA"),
    "PHA-DEX"=c("PHA", "PHA-DEX"))
contrast_nn <- names(contrast_ls)
metaNew <- map_dfc(contrast_nn, function(nn){
   ### 
   oneX <- contrast_ls[[nn]] 

   ### shrinkage
   res2 <- res%>%filter(MCls==oneMCl, contrast==nn, qval<0.1, abs(beta)>th0)
   s0 <- rep(0, length(rn))
   names(s0) <- rn
   s0[res2$gene] <- 1 
    
   ##  lfc <- abs(res2$beta)
   ## lfc[is.na(lfc)] <- 0  
   ## names(lfc) <- res2$gene
    
   ### matrix 1 
   ii <- meta%>%filter(treats==oneX[1])%>%dplyr::pull(NEW_BARCODE)
   X1 <-X[,ii]
   mu1 <- rowMeans(X1)

   ### matrix 2
   ii <- meta%>%filter(treats==oneX[2])%>%dplyr::pull(NEW_BARCODE)
   X2 <-X[,ii]
   mu2 <- rowMeans(X2)

   ## Difference
   Diff <- as.matrix((mu2-mu1)*(1/var)*s0)
   Diff[is.na(Diff)] <- 0
    
   ##
   z1 <- as.vector(crossprod(Xc, Diff))
   z1
 })
    
metaNew <- as.data.frame(metaNew) 
names(metaNew) <- paste("z_", contrast_nn, sep="")  
metaNew <- metaNew%>%
    mutate(NEW_BARCODE=colnames(X), treats=gsub(".*-ATAC-|_.*", "", NEW_BARCODE))
 
###
opfn <- paste(outdir, "1_MCls.", oneMCl, ".DLDA.rds", sep="")
write_rds(metaNew, opfn)
    
}###



#########################
### summary LDA plots ###
#########################

### scatter plots plus density plots
LDAplot <- function(df, labXY=c("x","y"), is.title=TRUE){
    
  coltreat <- c("CTRL"="#828282",
                "LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
                "PHA"="#a6cee3", "PHA-DEX"="#1f78b4")
    
  MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")

  fig_ls <- lapply(MCls, function(oneMCl){
    df2 <- df%>%filter(MCls==oneMCl)
    fig0 <- ggplot(df2, aes(x=x, y=y))+
      rasterise(geom_point(aes(colour=factor(treats)), size=0.1), dpi=300)+
      scale_colour_manual(values=coltreat, guide=guide_legend(override.aes=list(size=2)))+
    xlab(labXY[1])+ylab(labXY[2])+
    ggtitle(oneMCl)+
    theme_bw()+
    theme(legend.position="none",
          axis.title=element_text(size=9),
          axis.text=element_text(size=7))
    ###
    if (is.title){
       fig0 <- fig0+theme(plot.title=element_text(hjust=0.5, vjust=-2, size=12))
    }else{
       fig0 <- fig0+theme(plot.title=element_blank())
    }
    ###
    fig0 <- ggMarginal(fig0, groupColour=T, groupFill=F, size=3)
    fig0
  })

  lab2 <- c("CTRL"="CTRL", "LPS"="LPS", "LPS-DEX"="LPS+DEX",
            "PHA"="PHA", "PHA-DEX"="PHA+DEX")
  legend2 <- get_legend(
    ggplot(df%>%filter(MCls=="Bcell"), aes(x, y))+
    geom_point(aes(colour=factor(treats)), size=0.1)+
    scale_colour_manual(values=coltreat, guide=guide_legend(override.aes=list(size=2)),labels=lab2)+
    theme_bw()+
    theme(legend.title=element_blank(),
          legend.background=element_rect(colour=NA, fill=NA),
          legend.text=element_text(size=8),
          legend.key=element_rect(fill=NA),
          legend.key.size=grid::unit(1,"lines"))
   ) ##end for legend2

  figures <- list(fig_ls=fig_ls, fig_legend=legend2)

  figures

} ###    
    
    

### Read data
MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
df <- map_dfr(MCls, function(oneMCl){
   ##
   cat(oneMCl, "\n")
   fn <- paste("./5_dynamic.outs/", "1_MCls.", oneMCl, ".DLDA.rds", sep="")
   dd <- read_rds(fn)%>%mutate(MCls=oneMCl)
   dd
})


################
### Figure 1 ###
################


## ### LPS vs PHA
## df1 <- df%>%
##    dplyr::select(NEW_BARCODE,treats, MCls, z_LPS, z_PHA)%>%
##    dplyr::rename(x=z_LPS, y=z_PHA)
## figures <- LDAplot(df1, labXY=c("DLDA of LPS", "DLDA of PHA"), is.title=TRUE)
## fig1 <- figures$fig_ls
## ## legend2 <- figures$fig_legend

## ## figfn <- paste("./9_RNA.dynamic2_output/Filter2_DEG6571/", option, "/Figure1.X1_LPSandPHA.png", sep="")
## ## png(figfn, width=900, height=700, res=130)
## ## fig1 <- plot_grid(fig_ls[[1]], fig_ls[[2]],
## ##                   fig_ls[[3]], fig_ls[[4]],
## ##                   nrow=2, ncol=2, align="hv",axis="tb")

## ## print(plot_grid(fig1, legend2, rel_widths=c(4,1)))
## ## dev.off()


## ### LPS-DEX vs PHA-DEX
## df2 <- df%>%
##    dplyr::select(NEW_BARCODE, treats, MCls, "z_LPS-DEX", "z_PHA-DEX")%>%
##    dplyr::rename(x="z_LPS-DEX", y="z_PHA-DEX")
## figures <- LDAplot(df2, labXY=c("DLDA of LPS+DEX", "DLDA of PHA-DEX"), is.title=FALSE)
## fig2 <- figures$fig_ls
## legend2 <- figures$fig_legend

## ## figfn <- paste("./7_dynamic.outs/Figure1.X2_LPS-DEXandPHA-DEX.png", sep="")
## ## png(figfn, width=900, height=700, res=130)
## ## fig1 <- plot_grid(fig_ls[[1]], fig_ls[[2]],
## ##                   fig_ls[[3]], fig_ls[[4]],
## ##                   nrow=2, ncol=2, align="hv",axis="tb")

## ## print(plot_grid(fig1, legend2, rel_widths=c(4,1)))
## ## dev.off()

## figfn <- "./5_dynamic.outs/Figure1_combLDA_LPS_PHA.pdf"
## pdf(figfn, width=12, height=6)
## comb <- plot_grid(fig1[[1]], fig1[[2]], fig1[[3]], fig1[[4]],
##    fig2[[1]], fig2[[2]], fig2[[3]], fig2[[4]], nrow=2, ncol=4,
##    label_y=0.8,
##    labels=c("A", "B", "C", "D", "", "", "", ""),
##    label_fontface="plain", align="hv")
## print(plot_grid(comb, legend2, rel_widths=c(6,1)))
## dev.off()      



################
### Figure 2 ###
################

### LPS vs LPS-DEX
df1 <- df%>%
   dplyr::select(NEW_BARCODE, treats,  MCls, z_LPS, "z_LPS-DEX")%>%
   dplyr::rename(x="z_LPS-DEX", y=z_LPS)%>%
   filter(treats%in%c("CTRL", "LPS", "LPS-DEX"))
###
figures <- LDAplot(df1, labXY=c("DLDA of LPS+DEX", "DLDA of LPS"))
fig1 <- figures$fig_ls

 
### PHA vs PHA-DEX
df2 <- df%>%
   dplyr::select(NEW_BARCODE, treats, MCls, z_PHA, "z_PHA-DEX")%>%
   dplyr::rename(x="z_PHA-DEX", y="z_PHA")%>%
   filter(treats%in%c("CTRL", "PHA", "PHA-DEX"))
figures <- LDAplot(df2, labXY=c("DLDA of PHA+DEX", "DLDA of PHA"), is.title=F)
fig2 <- figures$fig_ls
legend2 <- figures$fig_legend


figfn <- "./5_dynamic.outs/Figure3.2_comb.LDA.pdf"
pdf(figfn, width=12, height=6)
comb <- plot_grid(fig1[[1]], fig1[[2]], fig1[[3]], fig1[[4]],
   fig2[[1]], fig2[[2]], fig2[[3]], fig2[[4]], nrow=2, ncol=4,
   label_y=0.8,
   labels=c("A", "B", "C", "D", "", "", "", ""),
   label_fontface="plain", align="hv")
print(plot_grid(comb, legend2, rel_widths=c(6,1)))
dev.off()




###############
### T-cells ###
###############

df <- read_rds("./5_dynamic.outs/1_MCls.Tcell.DLDA.rds")
coltreat <- c("CTRL"="#828282", "LPS"="#fb9a99", "LPS-DEX"="#e31a1c", "PHA"="#a6cee3", "PHA-DEX"="#1f78b4")

df2 <- df%>%
   dplyr::select(NEW_BARCODE, treats,  z_LPS, "z_LPS-DEX")%>%
   dplyr::rename(x="z_LPS-DEX", y="z_LPS")%>%
   filter(treats%in%c("CTRL", "LPS", "LPS-DEX"))

p1 <- ggplot(df2, aes(x=x, y=y))+
   rasterise(geom_point(aes(colour=factor(treats)), size=0.1), dpi=300)+
      scale_colour_manual(values=coltreat, guide=guide_legend(override.aes=list(size=2)))+
    xlab("DLDA of LPS+DEX")+ylab("DLDA of LPS")+
    theme_bw()+
    theme(legend.position="none",
          axis.title=element_text(size=14),
          axis.text=element_text(size=10))
p1 <- ggMarginal(p1, groupColour=T, groupFill=F, size=3)


###
###

df3 <- df%>%
   dplyr::select(NEW_BARCODE, treats,  z_PHA, "z_PHA-DEX")%>%
   dplyr::rename(x="z_PHA-DEX", y="z_PHA")%>%
   filter(treats%in%c("CTRL", "PHA", "PHA-DEX"))

p2 <- ggplot(df3, aes(x=x, y=y))+
   rasterise(geom_point(aes(colour=factor(treats)), size=0.1), dpi=300)+
      scale_colour_manual(values=coltreat, guide=guide_legend(override.aes=list(size=2)))+
    xlab("DLDA of PHA+DEX")+ylab("DLDA of PHA")+
    theme_bw()+
    theme(legend.position="none",
          axis.title=element_text(size=14),
          axis.text=element_text(size=10))
p2 <- ggMarginal(p2, groupColour=T, groupFill=F, size=3)

###
###
figfn <- "./5_dynamic.outs/Figure3.3_comb.DLDA.png"
png(figfn, width=420, height=800, res=120)
plot_grid(p1, p2, nrow=2, align="hv", axis="lr")
dev.off()


#############################################
### Dynamical changes of motif activities ###
#############################################

### sliding window
slideFun <- function(X_sort, lda, win=0.1, step=0.001){
###
  win <- trunc(ncol(X_sort)*win)
  step <- trunc(ncol(X_sort)*step)
    
  X <- X_sort  
  nlen <- ncol(X) 

  Xnew <- NULL
  ldaNew <- NULL  
  s0 <- 1
  while(TRUE){
    ##
    s1 <- s0+win-1
    if (s1>nlen) break
    xi <- apply(X[,s0:s1], 1, sum)
    Xnew <- cbind(Xnew, xi)
    ldai <- mean(lda[s0:s1])
    ldaNew <- c(ldaNew, ldai)  
    s0 <- s0+step
  }
  res <- list(Xnew=Xnew, ldaNew=ldaNew)
  res  
}

###
getDynamicalMat <- function(X, meta, contrast, win=0.1, step=0.001){

 ##  treat0 <- contrast[1]
##   treat1 <- contrast[2]
## ###
##   barcode <- colnames(X)
##   x1 <- X[,grepl(treat0, barcode)]
##   mu1 <- apply(x1, 1, mean)
##   x2 <- X[,grepl(treat1, barcode)]
##   mu2 <- apply(x2, 1, mean)
##   diff <- data.frame(gene=gsub("S-|\\..*", "", rownames(X)), beta=mu2-mu1)%>%
##     arrange(beta) 
  ## geneSel <- diff$gene
  ## mat2 <- mat2[geneSel,]

### re-order cell
  meta_sort <- meta%>%arrange(z)

## treatment 1 
  treat1 <- contrast[2]
  ##cat(treat1, "\n")  
  meta1 <- meta_sort%>%filter(treats==treat1)  
  mat1 <- as.matrix(X[, meta1$NEW_BARCODE])
  #rownames(mat2) <- gsub("S-|\\..*", "", rownames(mat2))
  tmp1 <- slideFun(mat1, lda=meta1$z, win=win, step=step)
  mat1 <- tmp1$Xnew
  lda1 <- tmp1$ldaNew  
### contrast treatment
  treat0 <- contrast[1]
  ##cat(treat0, "\n")  
  meta0 <- meta_sort%>%filter(treats==treat0)  
  mat0 <- as.matrix(X[, meta0$NEW_BARCODE])
  #rownames(mat2) <- gsub("S-|\\..*", "", rownames(mat2))
  tmp0 <- slideFun(mat0, lda=meta0$z, win=win, step=step)
  mat0 <- tmp0$Xnew
  lda0 <- tmp0$ldaNew  
###
  min1 <- apply(mat1, 1, min)
  max1 <- apply(mat1, 1, max)
  Range <- max1-min1
  mat1 <- sweep(mat1, 1, min1, "-")
  mat1 <- sweep(mat1, 1, Range, "/")  
    
  ## mMax <- apply(cbind(mat1, mat0), 1, max)
    
  ## mat1 <- mat1[mMax>0,]
  ## mat0 <- mat0[mMax>0,]  
  ## mMax <- mMax[mMax>0]  
  ## mat1 <- sweep(mat1, 1, mMax, "/")
  ## mat1 <- as.matrix(mat1)
  ## ##
  ## mat0 <- sweep(mat0, 1, mMax, "/")
  ## mat0 <- as.matrix(mat0)
###
  mat_ls <- list(mat1=mat1, lda1=lda1)  
  mat_ls  
###    
}


###setting colors
col_fun <-  colorRamp2(seq(0, 1, len=11), rev(brewer.pal(11, "Spectral")))


###
res <- read_rds("./2_motif.activities.outs/3_motif.diff.results.rds")
th0 <- quantile(abs(res$beta),probs=0.9)

##
## MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
## LDA <- c("LPS", "LPS-DEX", "PHA", "PHA-DEX")
## dataset <- data.frame(MCls=rep(MCls, each=4), LDA=rep(LDA, times=4))



contrast_ls <- list("LPS"=c("CTRL", "LPS"),
   "LPS-DEX"=c("LPS", "LPS-DEX"),
   "PHA"=c("CTRL", "PHA"),
   "PHA-DEX"=c("PHA", "PHA-DEX"))
coltreat <- c("CTRL"="#828282",
   "LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
   "PHA"="#a6cee3", "PHA-DEX"="#1f78b4")

oneMCl <- "Tcell"
lda <- "PHA-DEX"
contrast <- contrast_ls[[lda]]
treat0 <- contrast[1]
treat1 <- contrast[2] 
ii <- 4
col0 <- coltreat[treat1]


###
### subset data
fn <- "./2_motif.activities.outs/1_scATAC.motifActivities.rds"
sc <- read_rds(fn)
motif <- Motifs(sc)
motif.names <- unlist(motif@motif.names)
###
sc2 <- subset(sc, subset=MCls==oneMCl)
X <- sc2@assays$chromvar@data
rn <- rownames(X)
rownames(X) <- motif.names[rn]

### sort cells by LDA
fn <- "./5_dynamic.outs/1_MCls.Tcell.DLDA.rds"
meta <- read_rds(fn)
meta2 <- meta%>%
   dplyr::select(NEW_BARCODE, treats)%>%mutate(z=meta[,ii])
meta2 <- meta2%>%filter(treats%in%c("PHA", "PHA-DEX"))

res0 <- res%>%filter(MCls==oneMCl, contrast==lda, abs(beta)>th0, qval<0.1)
    
###
### heatmap for treat1 by cluster genes
Xi <- X[as.character(res0$motif),meta2$NEW_BARCODE]
xmin <- apply(Xi, 1, min)
X2 <- sweep(Xi, 1, xmin, "-")
sd2 <- apply(X2, 1, sd)
X2 <- sweep(X2, 1, sd2, "/")
##
mat_ls <- getDynamicalMat(X2, meta2, contrast, win=0.1, step=0.001)
mat1 <- mat_ls[[1]]
###
## rownames(mat1) <- motif.names[rownames(mat1)]

### heatmap
### column annotation
z <- mat_ls$lda1
z2 <- (z-min(z))/(max(z)-min(z))
#z3 <- z2[z2>0.23] 
#breaks <- c(0,quantile(z3,probs=seq(0,1,length.out=99))) 
col2 <- colorRamp2(seq(0,1,length.out=100),
   colorRampPalette(c("white", coltreat[[treat1]]))(100))

col_ha <- HeatmapAnnotation(pseudotime=z2, col=list(pseudotime=col2),
   show_legend=FALSE,                         
   annotation_name_gp=gpar(fontsize=8),                         
   simple_anno_size=unit(0.3,"cm") )

###
# <- colorRampPalette(rev(brewer.pal(n=7, name="RdBu")))(100) 
## b <- res0$beta
## th0 <- as.numeric(quantile(abs(b),probs=0.99)) 
## b2 <- b[abs(b)<th0] ### quantile(abs(b),probs=0.99)
## breaks <- c(min(b),quantile(b2,probs=seq(0,1,length.out=98)),max(b)) 
## col2 <- colorRamp2(breaks,
##    colorRampPalette(rev(brewer.pal(n=7, name="RdBu")))(100) ) 
## row_ha <- rowAnnotation(LFC=res0$beta, col=list(LFC=col2),
##    annotation_legend_param=list(LFC=list(grid_width=unit(0.3,"cm"),
##    labels_gp=gpar(fontsize=8), title_gp=gpar(fontsize=10), title="LFC")),    
##    width=unit(0.5,"cm"), annotation_name_gp=gpar(fontsize=8)) 


fig1 <- Heatmap(mat1, col=col_fun,
   cluster_rows=TRUE, cluster_columns=FALSE,
   row_km=2, show_parent_dend_line=FALSE,
   show_row_names=T,
   row_names_side="right",
   row_names_gp=gpar(fontsize=7),
   show_row_dend=T,
   show_column_names=FALSE,
   column_title=paste("pseudotime(", treat1, ")", sep=""),
   column_title_gp=gpar(fontsize=10),
   top_annotation=col_ha,
   heatmap_legend_param=list(title="Relative activities",
   title_gp=gpar(fontsize=8),
   labels_gp=gpar(fontsize=8), grid_width=unit(0.3, "cm")),
   use_raster=TRUE, raster_device="png")

 

figfn <- "./5_dynamic.outs/Figure4_Tcell_ldaPHA-DEX_trtPHA-DEX.png"
png(figfn, height=600, width=500, res=120)
set.seed(0)
fig1 <- draw(fig1)
r.list <- row_order(fig1)
r.dend <- row_dend(fig1) 
dev.off()


###cluster genes
clu_df <- lapply(names(r.list),function(i){
   out <- data.frame(genes=rownames(mat1)[r.list[[i]]],
                     Cluster=i, stringsAsFactors=FALSE)
   out
})%>%do.call(rbind,.)
rownames(clu_df) <- clu_df$gene
###
opfn <- paste(outdir, "6.3_Tcell_ldaPHA-DEX_trtPHA-DEX.cluster.csv", sep="")
write.csv(clu_df, opfn, row.names=FALSE)
