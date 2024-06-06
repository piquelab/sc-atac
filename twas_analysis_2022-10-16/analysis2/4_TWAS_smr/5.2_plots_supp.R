##
library(Matrix)
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(data.table)
library(qvalue)
library(annotables)

##
## library(GenomicRanges)
## library(Seurat)
## library(SeuratDisk)
## library(SeuratData)
## library(Signac)  ##, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
## library(SeuratWrappers)
## library(cicero, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
## library(monocle3)


##library(ComplexHeatmap)
library(circlize)
library(openxlsx)
library(cowplot)
library(ggrepel)
library(scales)
library(ReactomePA)
##

rm(list=ls())


###
###
outdir <- "./5_pub.outs/2_supp_plots/"
if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)



 
#############################################
### histogram distribution of p-value  
#############################################


trait <- sort(read.table("traits.txt")$V1)[3]

methods <- c("aloft_minP", "aloft_topPIP_union", "gtex_minP", "gtex_topPIP_union")


###
### plot data 
plotDF <- map_dfr(methods, function(ii){
    ###
    fn <- paste("./twas_smr.outs/", trait, "_", ii, "_twas.txt.gz", sep="")
    res <- fread(fn, header=T, data.table=F)
    ###
    pval <- res$pval_gwas
    fdr <- qvalue(pval) 
    ##
    res$FDR <- fdr$qvalues
    res$conditions <- ii
    ##
    res2 <- res%>%dplyr::select(gene, pval_gwas, conditions)
    res2
})
##


###
### annotation 
annoDF <- map_dfr(methods, function(ii){
   ###
   fn <- paste("./twas_smr.outs/", trait, "_", ii, "_twas.txt.gz", sep="")
   res <- fread(fn, header=T, data.table=F)
   ###
   pval <- res$pval_gwas
   fdr <- qvalue(pval)    

   res$FDR <- fdr$qvalues    
   res2 <- res%>%filter(FDR<0.1)

   ###
   pi0 <- round(fdr$pi0, digits=3)
   ngene0 <- length(unique(res$gene))
    
   eq2 <- as.expression(bquote(~pi==.(pi0)))
    
   cat(ii, length(unique(res2$gene)), "\n")
   df2 <- tibble(conditions=ii,
                 nsigs=length(unique(res2$gene)), pi=pi0, eq=eq2,
                 yline=ngene0*pi0, ngene=ngene0)
   ##    
   df2
})

##
##
annoDF <- annoDF%>%mutate(yline2=yline/40, xpos=0.75, ypos=1500)

facet_lab <- as_labeller(c("aloft_minP"="ALOFT_SMR_standard", "aloft_topPIP_union"="ALOFT_SMR_annotation",
              "gtex_minP"="GTEx_SMR_standard", "gtex_topPIP_union"="GTEx_SMR_annotation"))
  
p0 <- ggplot()+
   geom_histogram(data=plotDF, aes(x=pval_gwas, color=conditions), fill=NA, bins=40)+
   geom_text(data=annoDF, aes(x=xpos, y=ypos, label=eq), size=4, parse=T)+
   geom_hline(data=annoDF, aes(yintercept=yline2), linetype="dashed")+ 
   scale_color_manual(values=c("aloft_minP"="#f1b6da", "aloft_topPIP_union"="#d01c8b",
                               "gtex_minP"="#b8e186", "gtex_topPIP_union"="#4dac26"))+
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
         strip.text=element_text(size=12))
 
figfn <- paste(outdir, "FigS5_1_asthma_pval.hist.pdf", sep="")
## ggsave(figfn, p2, width=580, height=550, units="px", dpi=120)
pdf(figfn, width=7, height=7)
print(p0)
dev.off()























## ###
## ### pathway enrichment analysis  
## trait <- sort(read.table("traits.txt")$V1)[3]

## ###
## ### topPIP, annotation
## fn <- paste("./4_INTACT.outs/", trait, "_aloft_topPIP_union_combinfor.txt", sep="")
## res <- read.table(fn, header=T, sep="\t")
## sig <- res%>%filter(FDR<0.1)%>%pull(Gene)%>%unique()

## ## sig <- res%>%filter(FDR<0.1, peaking_d==3)%>%pull(Gene)%>%unique()


## ###
## ### minP, w/o annotation
## fn <- paste("./4_INTACT.outs/", trait, "_aloft_minP_intact.txt", sep="")
## res2 <- read.table(fn, header=T)
## sig2 <- res2%>%filter(FDR<0.1)%>%pull(Gene)%>%unique()


## ###
## ### background genes 
## BG <- unique(sig)
## geneBG <- bitr(BG, fromType="ENSEMBL", toType="ENTREZID", OrgDb=org.Hs.eg.db)


## ###
## ### shared genes
## olap <- intersect(sig, sig2)
## df0 <- data.frame(ENSEMBL=olap)%>%
##     mutate(cluster="CL1_share")%>%
##     inner_join(geneBG, by="ENSEMBL")


## ###
## ### unique 
## unq2 <- setdiff(sig, sig2)
## df2 <- data.frame(ENSEMBL=unq2)%>%
##     mutate(cluster="CL2_only")%>%
##     inner_join(geneBG, by="ENSEMBL")


## ## ###
## ## df3 <- data.frame(ENSEMBL=unique(sig))%>%
## ##     mutate(cluster="CL3_annot")%>%
## ##     inner_join(geneBG, by="ENSEMBL")

## geneCluster <- rbind(df0, df2) 

## ###
## ### GO enrichment 

## cluster_gr <- sort(unique(geneCluster$cluster))
## res <- map_dfr(cluster_gr, function(ii){
## ###    
## ###    
##    df0 <- geneCluster%>%filter(cluster==ii) 
##    res0 <- enrichGO(gene=df0$ENTREZID, universe=geneBG$ENTREZID,
##                OrgDb="org.Hs.eg.db", ont="ALL", pvalueCutoff=1, qvalueCutoff=1, minGSSize=0, maxGSSize=1000)
              
##    res0 <- res0%>%mutate(maxGSSize=as.numeric(gsub("/.*", "", BgRatio)), padj2=p.adjust(pvalue, method="BH"))%>%
##        as.data.frame()
##    res0$cluster <- ii 
##    cat(ii, sum(res0$padj2<0.1), "\n") 
##    res0 
## })

## opfn <- paste(outdir, "1_GO_enriched.rds", sep="")
## write_rds(res, file=opfn)


## ###
## ### KEGG enrichment
## cluster_gr <- sort(unique(geneCluster$cluster))
## res <- map_dfr(cluster_gr, function(ii){
## ###    
## ###    
##    df0 <- geneCluster%>%filter(cluster==ii) 
##    res0 <- enrichKEGG(gene=df0$ENTREZID, universe=geneBG$ENTREZID,
##         pvalueCutoff=1, qvalueCutoff=1, minGSSize=0, maxGSSize=1000)
##    ## 
##    res0 <- res0%>%mutate(maxGSSize=as.numeric(gsub("/.*", "", BgRatio)),
##                          padj2=p.adjust(pvalue, method="BH"))%>%
##        as.data.frame()
##    res0$cluster <- ii 
    
##    cat(ii, sum(res0$padj2<0.1), "\n") 
##    res0 
## })

## opfn <- paste(outdir, "2_KEGG_enriched.rds", sep="")
## write_rds(res, file=opfn)    


## ###
## ### Reactome enrichment analysis
## cluster_gr <- sort(unique(geneCluster$cluster))
## res <- map_dfr(cluster_gr, function(ii){
## ###    
## ###    
##    df0 <- geneCluster%>%filter(cluster==ii) 
##    res0 <- enrichPathway(gene=df0$ENTREZID, universe=geneBG$ENTREZID,
##         pvalueCutoff=1, qvalueCutoff=1, minGSSize=0, maxGSSize=1000)
##    ## 
##    res0 <- res0%>%mutate(maxGSSize=as.numeric(gsub("/.*", "", BgRatio)),
##                          padj2=p.adjust(pvalue, method="BH"))%>%
##        as.data.frame()
##    res0$cluster <- ii 
    
##    cat(ii, sum(res0$padj2<0.1), "\n") 
##    res0 
## })

## opfn <- paste(outdir, "3_Reactome_enriched.rds", sep="")
## write_rds(res, file=opfn)    



## ##########################################
## ### extract significant results 
## ###########################################

## ### GO
## fn <- paste(outdir, "1_GO_enriched.rds", sep="")
## res <- read_rds(fn)
## res_sig <- res%>%
##     filter(pvalue<0.05, maxGSSize>=5, maxGSSize<500)%>%group_by(cluster)%>%
##     arrange(pvalue, .by_group=T)
## ## res_sig%>%group_by(cluster)%>%summarize(nsig=n(), .groups="drop")

## opfn <- paste(outdir, "1.2_sig_GO_enriched.xlsx", sep="")
## write.xlsx(res_sig, file=opfn)

## ### KEGG
## fn <- paste(outdir, "2_KEGG_enriched.rds", sep="")
## res <- read_rds(fn)
## res_sig <- res%>%
##     filter(pvalue<0.05, maxGSSize>=5, maxGSSize<500)%>%group_by(cluster)%>%
##     arrange(pvalue,.by_group=T)
## ## res_sig%>%group_by(cluster)%>%summarize(nsig=n(), .groups="drop")

## opfn <- paste(outdir, "2.2_sig_KEGG_enriched.xlsx", sep="")
## write.xlsx(res_sig, file=opfn)

 
## ###
## fn <- paste(outdir, "3_Reactome_enriched.rds", sep="")
## res <- read_rds(fn)
## res_sig <- res%>%
##     filter(pvalue<0.05, maxGSSize>=5, maxGSSize<500)%>%group_by(cluster)%>%
##     arrange(pvalue, .by_group=T)
## res_sig%>%group_by(cluster)%>%summarize(nsig=n(), .groups="drop")

## opfn <- paste(outdir, "3.2_sig_Reactome_enriched.xlsx", sep="")
## write.xlsx(res_sig, file=opfn)

## ###
## ### enrichment analysis
## ## geneCluster <- rbind(df0, df2, df3)
## ## cg <- compareCluster(ENTREZID~cluster, data=geneCluster, universe=geneBG$ENTREZID,
## ##    fun="enrichGO", OrgDb="org.Hs.eg.db",
## ##    pvalueCutoff=1, qvalueCutoff=1, ont="ALL", minGSSize=0, maxGSSize=1000)                  

## ## ###
## ## opfn <- paste(outdir, "1_asthma_risk_genes_enrichGO.rds", sep="")
## ## write_rds(cg, file=opfn)
 


## #####################################
## ### GO enrichment analysis 
## #####################################

## ## outdir <- "./5_pub.outs/2_supp_plots/"
## ## ##
## ## fn <- paste(outdir, "1_asthma_risk_genes_enrichGO.rds", sep="")
## ## cg <- read_rds(fn)


## for ( ii in  sort(unique(geneCluster$cluster))){
## ###    
## ###
##    cat(ii, "\n")
    
##    df0 <- geneCluster%>%filter(cluster==ii) 
##    cg <- enrichGO(gene=df0$ENTREZID, universe=geneBG$ENTREZID,
##                OrgDb="org.Hs.eg.db", ont="ALL", pvalueCutoff=1, qvalueCutoff=1, minGSSize=0, maxGSSize=1000)
              
##    cg <- cg%>%mutate(maxGSSize=as.numeric(gsub("/.*", "", BgRatio)), padj2=p.adjust(pvalue, method="BH"))
  
##    cg2 <- cg%>%filter(maxGSSize>=5, maxGSSize<500)%>%arrange(pvalue)%>%as.data.frame()
        
##    plotDF <- cg2[1:10,]%>%
##      dplyr::select(Description, Count, pvalue, padj2)%>%
##      mutate(log10p=-log10(pvalue),
##             Descrip2=fct_reorder(Description, Count))
    
##      cat(ii, sum(cg2$padj2<0.1), "\n") 
    
## ###
## ### setting colors
## ## y0 <- plotDF$log10p
## ## mybreaks <- quantile(y0, probs=seq(0, 1, length.out=10))
## ## mycol <- colorRampPalette(c("blue", "red"))(10)

## ## mycol <- colorRamp2(mybreaks, colorRampPalette(c("blue", "red"))(10)) 

 
## p0 <- ggplot(plotDF, aes(x=Descrip2, y=Count, fill=log10p))+
##     geom_bar(stat="identity")+
##     scale_fill_gradient(name=bquote(-log[10]~italic(p)),
##         low="blue", high="red", n.breaks=5,
##         ## colours=mycol, values=seq(0, 1, length.out=10),
##         guide=guide_colourbar(barwidth=grid::unit(0.4, "cm"),
##                               barheight=grid::unit(4, "cm")))+
##     scale_x_discrete(labels=label_wrap(55))+
##     ylab("#Genes in GO terms ")+
##     ggtitle(ii)+
##     coord_flip()+
##     theme_bw()+
##     theme(plot.title=element_text(hjust=0.5, size=10),
##           axis.title.x=element_text(size=10),
##           axis.title.y=element_blank(),
##           axis.text.x=element_text(size=8),
##           axis.text.y=element_text(size=8),
##           legend.title=element_text(size=8),
##           legend.text=element_text(size=8))

## figfn <- paste(outdir, "Figure1_GO_", ii, ".bar.png", sep="")
## ggsave(figfn, p0, width=680, height=420, units="px", dpi=120)

## }


## ##################
## ### KEGG 
## ##################

## for ( ii in  sort(unique(geneCluster$cluster))){
## ###    
## ###
##    cat(ii, "\n")
    
##    df0 <- geneCluster%>%filter(cluster==ii) 
##    cg <- enrichKEGG(gene=df0$ENTREZID, universe=geneBG$ENTREZID,
##         pvalueCutoff=1, qvalueCutoff=1, minGSSize=0, maxGSSize=1000)
              
##    cg <- cg%>%mutate(maxGSSize=as.numeric(gsub("/.*", "", BgRatio)), padj2=p.adjust(pvalue, method="BH"))
  
##    cg2 <- cg%>%filter(maxGSSize>=5, maxGSSize<500)%>%arrange(pvalue)%>%as.data.frame()
        
##    plotDF <- cg2[1:10,]%>%
##      dplyr::select(Description, Count, pvalue, padj2)%>%
##      mutate(log10p=-log10(pvalue),
##             Descrip2=fct_reorder(Description, Count))
    
##   cat(ii, sum(cg2$padj2<0.1), "\n") 
## ###
## ### setting colors
## ## y0 <- plotDF$log10p
## ## mybreaks <- quantile(y0, probs=seq(0, 1, length.out=10))
## ## mycol <- colorRampPalette(c("blue", "red"))(10)

## ## mycol <- colorRamp2(mybreaks, colorRampPalette(c("blue", "red"))(10)) 

 
## p0 <- ggplot(plotDF, aes(x=Descrip2, y=Count, fill=log10p))+
##     geom_bar(stat="identity")+
##     scale_fill_gradient(name=bquote(-log[10]~italic(p)),
##         low="blue", high="red", n.breaks=5,
##         ## colours=mycol, values=seq(0, 1, length.out=10),
##         guide=guide_colourbar(barwidth=grid::unit(0.4, "cm"),
##                               barheight=grid::unit(4, "cm")))+
##     scale_x_discrete(labels=label_wrap(55))+
##     ylab("#Genes in GO terms ")+
##     ggtitle(ii)+
##     coord_flip()+
##     theme_bw()+
##     theme(plot.title=element_text(hjust=0.5, size=10),
##           axis.title.x=element_text(size=10),
##           axis.title.y=element_blank(),
##           axis.text.x=element_text(size=8),
##           axis.text.y=element_text(size=8),
##           legend.title=element_text(size=8),
##           legend.text=element_text(size=8))

## figfn <- paste(outdir, "Figure2_KEGG_", ii, ".bar.png", sep="")
## ggsave(figfn, p0, width=680, height=420, units="px", dpi=120)

## }


## ##########################################
## ### Reactome enrichment analysis 
## ###########################################


 
## for ( ii in  sort(unique(geneCluster$cluster))){
## ###    
## ###
    
##    df0 <- geneCluster%>%filter(cluster==ii) 
##    cg <- enrichPathway(gene=df0$ENTREZID, universe=geneBG$ENTREZID,
##         pvalueCutoff=1, qvalueCutoff=1, minGSSize=0, maxGSSize=1000)
              
##    cg <- cg%>%mutate(maxGSSize=as.numeric(gsub("/.*", "", BgRatio)), padj2=p.adjust(pvalue, method="BH"))
  
##    cg2 <- cg%>%filter(maxGSSize>=5, maxGSSize<500)%>%arrange(pvalue)%>%as.data.frame()
        
##    plotDF <- cg2[1:10,]%>%
##      dplyr::select(Description, Count, pvalue, padj2)%>%
##      mutate(log10p=-log10(pvalue),
##             Descrip2=fct_reorder(Description, Count))

##    cat(ii, sum(cg2$padj2<0.1), "\n") 
## ###
## ### setting colors
## ## y0 <- plotDF$log10p
## ## mybreaks <- quantile(y0, probs=seq(0, 1, length.out=10))
## ## mycol <- colorRampPalette(c("blue", "red"))(10)

## ## mycol <- colorRamp2(mybreaks, colorRampPalette(c("blue", "red"))(10)) 

 
## p0 <- ggplot(plotDF, aes(x=Descrip2, y=Count, fill=log10p))+
##     geom_bar(stat="identity")+
##     scale_fill_gradient(name=bquote(-log[10]~italic(p)),
##         low="blue", high="red", n.breaks=5,
##         ## colours=mycol, values=seq(0, 1, length.out=10),
##         guide=guide_colourbar(barwidth=grid::unit(0.4, "cm"),
##                               barheight=grid::unit(4, "cm")))+
##     scale_x_discrete(labels=label_wrap(55))+
##     ylab("#Genes in GO terms ")+
##     ggtitle(ii)+
##     coord_flip()+
##     theme_bw()+
##     theme(plot.title=element_text(hjust=0.5, size=10),
##           axis.title.x=element_text(size=10),
##           axis.title.y=element_blank(),
##           axis.text.x=element_text(size=8),
##           axis.text.y=element_text(size=8),
##           legend.title=element_text(size=8),
##           legend.text=element_text(size=8))

## figfn <- paste(outdir, "Figure3_Reactome_", ii, ".bar.png", sep="")
## ggsave(figfn, p0, width=680, height=420, units="px", dpi=120)

## }

