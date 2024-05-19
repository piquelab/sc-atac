###
###
library(Matrix)
library(tidyverse)
library(data.table)
##library(clusterProfiler)
##library(org.Hs.eg.db)

## library(qvalue)
## library(annotables)

## ##
## library(GenomicRanges)
## library(Seurat)
## library(SeuratDisk)
## library(SeuratData)
## library(Signac)  ##, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
## library(SeuratWrappers)
## library(cicero, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
## library(monocle3)


##library(ComplexHeatmap)
##library(circlize)
## library(openxlsx)
library(cowplot)
library(ggrepel)
library(ggrastr)
library(grid)
library(lattice)
library(ggplotify)
##

rm(list=ls())

###
###
outdir <- "./3_plots.outs/"
if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)

### passing argument
args=commandArgs(trailingOnly=T)
if ( length(args)>0){
   ##
   trait <- args[1]
}else{
   trait <- "Asthma_Bothsex_afr_inv_var_meta_GBMI_052021_nbbkgt1"
}   


traits <- read.table("traits.txt")$V1

#########################
### read twas results ###
#########################


if (grepl("GBMI", trait)){
    ##
    fn <- paste("./2_impute.outs/", trait, "_impute.txt.gz", sep="")
    res <- fread(fn, header=T, data.table=F)
    res <- res[,c(1,2,6)]
    names(res) <- c("chr", "pos", "pval")
    trait0 <- gsub("_inv_.*", "", trait)
}else{
    ##
    fn <- paste("./2_impute.outs/", trait, "_impute.txt.gz", sep="")
    res <- fread(fn, header=T, data.table=F)
    res <- res[,c(1,2,6)]
    names(res) <- c("chr", "pos", "pval")
    res <- res%>%mutate(chr=as.integer(chr), pos=as.integer(pos))
    trait0 <- "Asthma_UKB"
}    

res$chr_pos <- paste(res$chr, res$pos, sep="_")


###
###
ngene <- nrow(res)
nchr <- max(res$chr)

###
res$BPcum <- NA
s <- 0
for (i in 1:nchr){
   x <- res[res$chr==i, "pos"]+s
   res[res$chr==i,"BPcum"] <- x
   s <- max(x)
}

axis.set <- res%>%
   group_by(chr)%>%
   summarize(midpoint=(max(BPcum)+min(BPcum))/2)

res <- res%>%mutate(log10p=-log10(pval), gr2=factor(chr%%2))

sig0 <- -log10(0.05/ngene)
sig1 <- -log10(0.01/ngene)




ymax <- max(res$log10p)+8

########################
### Manhattan plots
##########################
p1 <- ggplot(res)+
   rasterize(geom_point(aes(x=BPcum, y=log10p, color=factor(gr2)), size=0.01), dpi=50)+
   geom_hline(yintercept=sig1, color="red", linetype="dashed")+ 
   scale_x_continuous("chromosome", label=axis.set$chr, breaks=axis.set$midpoint,
      limits=c(min(res$BPcum), max(res$BPcum)) )+
   scale_y_continuous(bquote(-log[10]~"("~italic(p)~")"), limits=c(0,ymax))+
   scale_color_manual(values=c("0"="#2171b5", "1"="#6baed6"))+
   ggtitle(trait0)+ 
   theme_bw()+
   theme(legend.position="none",
         axis.title=element_text(size=10),
         axis.text.x=element_text(size=10, angle=90, vjust=0.5, hjust=1),
         axis.text.y=element_text(size=10),
         plot.title=element_text(hjust=0.5, size=12),
         panel.grid.major=element_blank(),
         panel.grid.minor=element_blank())
## p1 <- as.grob(p1)




###############
### qq plots
################
ci <- 0.95 
plotDF <- res%>%arrange(pval)%>%
    mutate(observed=-log10(pval), expected=-log10(ppoints(ngene)),
           clower=-log10(qbeta(p=(1-ci)/2, shape1=seq(ngene), shape2=rev(seq(ngene)))),
           cupper=-log10(qbeta(p=(1+ci)/2, shape1=seq(ngene), shape2=rev(seq(ngene)))) )
           
p2  <- ggplot(plotDF, aes(x=expected, y=observed))+
   rasterize(geom_point(size=0.01), dpi=50)+
   ## geom_ribbon(aes(ymax=cupper, ymin=clower), fill="grey", alpha=0.5)+ 
   geom_abline(slope=1, intercept=0)+
   xlab(bquote("Expected"~ -log[10]~"("~italic(plain(P))~")"))+
   ylab(bquote("Observed"~-log[10]~"("~italic(plain(P))~")"))+
   ggtitle(trait0)+ 
   theme_bw()+
   theme(legend.position="none",
         axis.text=element_text(size=10),
         axis.title=element_text(size=10),
         plot.title=element_text(hjust=0.5, size=12),
         panel.grid.major=element_blank(),
         panel.grid.minor=element_blank())

## p2 <- as.grob(p2)


plot_comb <- plot_grid(p1, p2, rel_widths=c(1.6, 1), nrow=1, ncol=2, aligh="h", axis="tb")
 

figfn <- paste(outdir, "Figure1_", trait, "_summary_comb.png", sep="")
ggsave(figfn, plot_comb, width=950, height=300, units="px", dpi=100)


###
### END

