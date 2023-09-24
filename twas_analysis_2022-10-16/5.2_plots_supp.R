###
###
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
outdir <- "./5_pub.outs/2_supp_plots/"
if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)



######################
### autosome genes ###
######################

autosome <- as.character(1:22) 
grch38_unq <- grch38%>%
    dplyr::filter(chr%in%autosome, grepl("protein", biotype))%>%
    distinct(ensgene, chr, .keep_all=T)%>%dplyr::select(gene=ensgene, chr, biotype, symbol)


#########################
### read twas results ###
#########################

fn <- "./Asthma_twas/ALOFT_topPIP_twas.txt.gz"
res <- fread(fn, header=T, data.table=F)
    
res <- res%>%dplyr::filter(gene%in%grch38_unq$gene)    
pval <- res$pval_gwas
fdr <- qvalue(pval)
res$FDR <- fdr$qvalues
 
res <- res%>%arrange(pval_gwas)%>%left_join(grch38_unq, by="gene")
res2 <- res%>%dplyr::filter(FDR<0.1)%>%dplyr::select(-gene_SNP, -biotype)


###
### torus annotation file
fn <- "../genetic.analysis_torus_2021-10-14/SNPannotation/4_SNPAnnot.outs/pct_0.02_combineNew/Union_torus.annot.gz"
anno_torus <- fread(fn, header=T, data.table=F)

## add annotation
res2 <- res2%>%left_join(anno_torus, by=c("genetic_variant"="SNP"))
 

###
### 
enloc <- read.table("./enloc_analysis/ALOFT_intact.txt", header=T)
enloc <- enloc%>%dplyr::select(gene=Gene, GLCP, PCG, FDR_intact=FDR)
res3 <- res2%>%left_join(enloc, by="gene")

###output
res3 <- res3%>%
    dplyr::select(gene, symbol, genetic_variant, chr, peaking_d, pval_gwas, FDR_twas=FDR,
    PIP, pval_eqtl, GLCP, PCG, FDR_intact) 
 
opfn <- gzfile(paste(outdir, "TableS5_1_asthma-risk-genes_ALOFT.txt.gz", sep=""))
write.table(res3, file=opfn, row.names=F, col.names=T, quote=F)

## opfn <- paste(outdir,  "TableS5_1_asthma-related-genes_ALOFT.xlsx", sep="")
## write.xlsx(res3, file=opfn, overwrite=T)

fn <- "./5_pub.outs/2_supp_plots/TableS5_1_asthma-risk-genes_ALOFT.txt.gz"
x <- fread(fn, header=T, data.table=F)



#########################################
### histogram distribution of p-value ###
#########################################

###
fn_ls <- list.files ("./Asthma_twas/", pattern="ALOFT|Whole") ## |Whole
fn_ls <- fn_ls[-2]

plotDF <- map_dfr(fn_ls, function(ii){
    ###
    fn <- paste("./Asthma_twas/", ii, sep="")
    res <- fread(fn, header=T, data.table=F)
    ##
    res <- res%>%dplyr::filter(gene%in%grch38_unq$gene)
    ###
    pval <- res$pval_gwas
    fdr <- qvalue(pval)
    ##
    res$FDR <- fdr$qvalues
    ##res2 <- res%>%filter(FDR<0.1)
    ##
    ii2 <- gsub("_twas.txt.gz", "", ii)
    res$conditions <- ii2
    ##
    res2 <- res%>%dplyr::select(gene, pval_gwas, conditions)
    res2
})
##


###
### annotation 
annoDF <- map_dfr(fn_ls, function(ii){
   ### 
   fn <- paste("./Asthma_twas/", ii, sep="")
   res <- fread(fn, header=T, data.table=F)
   res <- res%>%filter(gene%in%grch38_unq$gene)

   pval <- res$pval_gwas
   fdr <- qvalue(pval)
   res$FDR <- fdr$qvalues
    
   res2 <- res%>%filter(FDR<0.1)
   ##
   ii2 <- gsub("_twas.txt.gz", "", ii)

   ###
   pi0 <- round(fdr$pi0, digits=3)
   ngene0 <- length(unique(res$gene))
    
   eq2 <- as.expression(bquote(~pi==.(pi0)))
    
   cat(ii, length(unique(res2$gene)), "\n")
   df2 <- tibble(conditions=ii2,
                 nsigs=length(unique(res2$gene)), pi=pi0, eq=eq2, xpos=0.75, ypos=780,
                 yline=ngene0*pi0, ngene=ngene0)
   ##    
   df2
})


## p <- ggplot(plotDF, aes(x=pval_gwas, color=conditions))+
##    geom_density()+
##    scale_color_manual(values=c("ALOFT"="#e77f79", "ALOFT_PIP"="#d73027",
##                                "Whole_Blood"="#91bfdb", "Whole_Blood_PIP"="#4575b4"))+
##    ## scale_fill_manual(values=c("ALOFT"="#e77f79", "ALOFT_PIP"="#d73027",
##    ##                             "Whole_Blood"="#91bfdb", "Whole_Blood_PIP"="#4575b4"))+
##    ##xlab(bquote(-log[10]~"("~italic(plain(P))~")"))+ 
##    xlab(bquote(~italic(plain(P))~"from gwas"))+
##    ylab("Density")+
##    theme_bw()+
##    theme(legend.title=element_blank(),
##          legend.text=element_text(size=8),
##          legend.key.size=grid::unit(0.8, "lines"),
##          axis.text=element_text(size=10))

## figfn <- paste(outdir, "Figure2.1_pval_density.png", sep="")
## png(figfn, width=580, height=420, res=120)
## print(p)
## dev.off()

##
##
annoDF <- annoDF%>%mutate(yline2=yline/40)

facet_lab <- as_labeller(c("ALOFT_minP"="ALOFT_WBL_minP", "ALOFT_topPIP"="ALOFT_WBL_topPIP",
              "Whole_Blood_minP"="GTEx_WBL_minP",
              "Whole_Blood_topPIP"="GTEx_WBL_topPIP"))
  
p2 <- ggplot()+
   geom_histogram(data=plotDF, aes(x=pval_gwas, color=conditions), fill=NA, bins=40)+
   geom_text(data=annoDF, aes(x=xpos, y=ypos, label=eq), size=5, parse=T)+
   geom_hline(data=annoDF, aes(yintercept=yline2), linetype="dashed")+ 
   scale_color_manual(values=c("ALOFT_minP"="#f1b6da", "ALOFT_topPIP"="#d01c8b",
                               "Whole_Blood_minP"="#b8e186", "Whole_Blood_topPIP"="#4dac26"))+
   facet_wrap(~conditions,
              labeller=facet_lab,
              nrow=2, ncol=2, scales="fixed")+ 
   xlab(bquote(~italic(plain(P))~"of TWAS"))+
   ylab("Number of genes")+
   theme_bw()+
   theme(## legend.title=element_blank(),
         ## legend.text=element_text(size=10),
         ## legend.key.size=grid::unit(0.8, "lines"),
         legend.position="none",
         axis.text=element_text(size=12),
         axis.title=element_text(size=12),
         strip.text=element_text(size=14))
 
figfn <- paste(outdir, "FigS5_1_pval_hist.pdf", sep="")
pdf(figfn, width=7, height=7)
print(p2)
dev.off()



#####
###



