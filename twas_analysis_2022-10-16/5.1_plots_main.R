##
##
library(Matrix)
library(tidyverse)
##library(clusterProfiler)
##library(org.Hs.eg.db)
library(data.table)
library(qvalue)
library(annotables)

##
library(GenomicRanges)
library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(Signac)  ##, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
library(SeuratWrappers)
## library(cicero, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
## library(monocle3)


##library(ComplexHeatmap)
##library(circlize)
library(openxlsx)
library(cowplot)
library(ggrepel)
##

rm(list=ls())

###
###
outdir <- "./5_pub.outs/1_main_plots/"
if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)




########################
### manhattan plots ####
########################
  

###
## fn_ls <- list.files ("./Asthma_twas/", pattern="ALOFT|Whole") ## |Whole
## ii <- fn_ls[3]

###
autosome <- as.character(1:22) 
grch38_unq <- grch38%>%
    dplyr::filter(chr%in%autosome, grepl("protein", biotype))%>%
    distinct(ensgene, chr, .keep_all=T)%>%dplyr::select(gene=ensgene, chr, start, biotype, symbol)


###
### read twas data 
fn <- "./Asthma_twas/ALOFT_topPIP_twas.txt.gz"
resRaw <- fread(fn, header=T, data.table=F)
res <- resRaw%>%filter(gene%in%grch38_unq$gene)
res <- res%>%mutate(FDR=qvalue(pval_gwas)$qvalues, is_sig=ifelse(FDR<0.1, 1, 0))%>%
   left_join(grch38_unq, by="gene")%>%mutate(chr=as.integer(chr))


### torus annotation file
fn <- "../genetic.analysis_torus_2021-10-14/SNPannotation/4_SNPAnnot.outs/pct_0.02_combineNew/Union_torus.annot.gz"
anno_torus <- fread(fn, header=T, data.table=F)

res <- res%>%left_join(anno_torus, by=c("genetic_variant"="SNP"))

###
### colocalization
enloc <- read.table("./enloc_analysis/ALOFT_intact.txt", header=T)
enloc <- enloc%>%dplyr::select(gene=Gene, GLCP, PCG, FDR_intact=FDR)


res <- res%>%left_join(enloc, by="gene")


###
###
## x <- res%>%arrange(pval_gwas)%>%dplyr::filter(FDR<0.1)



###
ngene <- nrow(res)
nchr <- max(res$chr)

###
res$BPcum <- NA
s <- 0
for (i in 1:nchr){
   x <- res[res$chr==i, "start"]+s
   res[res$chr==i,"BPcum"] <- x
   s <- max(x)
}

axis.set <- res%>%
   group_by(chr)%>%
   summarize(midpoint=(max(BPcum)+min(BPcum))/2)

res <- res%>%mutate(log10p=-log10(pval_gwas), gr2=factor(chr%%2))


sig <- res%>%filter(is_sig==1)%>%pull(log10p)%>%min()

## x <- res%>%dplyr::filter(FDR<0.1)%>%arrange(pval_gwas)
###
### Manhattan plots

x <- res%>%dplyr::filter(is_sig==1)%>%arrange(pval_gwas)
x1 <- x[1:20,]
x2 <- x%>%dplyr::filter(FDR_intact<0.1, peaking_d==3)
annoDF2 <- rbind(x1,x2)%>%distinct(symbol, .keep_all=T)


ymax <- max(res$log10p)+8

p2 <- ggplot(res) +
   geom_point(aes(x=BPcum, y=log10p, color=factor(gr2), size=log10p))+
   geom_text_repel(data=annoDF2,
       aes(x=BPcum, y=log10p, label=symbol), box.padding=0.2,
            max.overlaps=15, fontface="italic", size=3)+      
   geom_hline(yintercept=sig, color="red", linetype="dashed")+ 
   scale_x_continuous("chromosome", label=axis.set$chr, breaks=axis.set$midpoint,
      limits=c(min(res$BPcum), max(res$BPcum)) )+
   scale_y_continuous(bquote(-log[10]~"("~italic(p)~")"), limits=c(0,ymax))+
   scale_color_manual(values=c("0"="#2171b5", "1"="#6baed6"))+
   scale_size_continuous(range=c(0.1,1.5)) +
   theme_bw()+
   theme(legend.position="none",
         axis.title=element_text(size=12),
         axis.text.x=element_text(size=10, angle=90, vjust=0.5, hjust=1),
         axis.text.y=element_text(size=12),
         panel.grid.major=element_blank(),
         panel.grid.minor=element_blank())

###
figfn <- paste(outdir, "Fig5_2_topPIP_asthma_manhattan.png", sep="")
png(figfn, width=900, height=450, pointsize=12, res=120)
print(p2)
dev.off()


geneSel <- c("GSDMB", "ORMDL3", "SLC25A46", "ERBB2", "BACH2",
             "IL4", "IKZF3", "FADS2", "WDR36", "ERBB3", "FLOT1", "MICB", "DEXI", "TRIM27", "TRIM26")

x2 <- x%>%filter(symbol%in%geneSel)
