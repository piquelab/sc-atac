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
    distinct(ensgene, chr, .keep_all=T)%>%dplyr::select(Gene=ensgene, start)


###
### read twas data 
fn <- "./5_pub.outs/2_supp_tables/TableS6_1_asthma-risk-genes_gtex.txt.gz"
res <- read.table(fn, header=T, sep="\t")
res <- res%>%left_join(grch38_unq, by="Gene")
names(res)[4] <- "chr"

res <- res%>%mutate(log10p=-log10(pval_twas), gr2=factor(chr%%2))

### 
### add x-axis information 
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
   summarize(midpoint=(max(BPcum)+min(BPcum))/2)%>%ungroup()



sig <- res%>%filter(FDR_twas<0.1)%>%pull(log10p)%>%min()


###
### annotation gene name
x <- res%>%dplyr::filter(FDR_twas<0.05)%>%arrange(pval_twas)
x1 <- x[1:20,]
x2 <- x%>%dplyr::filter(FDR_intact<0.05, peaking_d==3)
x2 <- x2[1:50,]
annoDF2 <- rbind(x1,x2)%>%distinct(symbol, .keep_all=T)


ymax <- max(res$log10p)+8
p1 <- ggplot(res)+
   geom_point(aes(x=BPcum, y=log10p, color=factor(gr2), size=log10p))+
   geom_text_repel(data=annoDF2,
       aes(x=BPcum, y=log10p, label=symbol), box.padding=0.2,
           min.segment.length=0, segment.curvature=1e-20, 
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
### QQ plots 
ci <- 0.95 
plotDF <- res%>%dplyr::select(Gene, pval=pval_twas)%>%arrange(pval)%>%
    mutate(observed=-log10(pval), expected=-log10(ppoints(ngene)),
           clower=-log10(qbeta(p=(1-ci)/2, shape1=seq(ngene), shape2=rev(seq(ngene)))),
           cupper=-log10(qbeta(p=(1+ci)/2, shape1=seq(ngene), shape2=rev(seq(ngene)))) )
           
p2  <- ggplot(plotDF, aes(x=expected, y=observed))+
   geom_point(color="#d01c8b", size=1)+
   ## geom_ribbon(aes(ymax=cupper, ymin=clower), fill="grey", alpha=0.5)+ 
   geom_abline(slope=1, intercept=0)+
   xlab(bquote("Expected"~ -log[10]~"("~italic(plain(P))~")"))+
   ylab(bquote("Observed"~-log[10]~"("~italic(plain(P))~")"))+
   theme_bw()+
   theme(legend.position="none",
         axis.text=element_text(size=12),
         axis.title=element_text(size=12),
         panel.grid.major=element_blank(),
         panel.grid.minor=element_blank())

## p2 <- as.grob(p2)


plot_comb <- plot_grid(p1, p2, rel_widths=c(1.7, 1),
    labels=c("A", "B"), label_fontface="plain", nrow=1, ncol=2, aligh="h", axis="tb")


###
figfn <- paste(outdir, "FigS5_2_asthma_manhattan_qq_comb.GTEx.png", sep="")
png(figfn, width=1200, height=450, pointsize=12, res=120)
print(plot_comb)
dev.off()

###
figfn <- paste(outdir, "FigS5_2_asthma_manhattan_qq_comb.GTEx.pdf", sep="")
pdf(figfn, width=12, height=4.5)
print(plot_comb)
dev.off()





###########################################################
### FigS5_3 Forest plots if enrichment in DEGs 
###########################################################




fn <- "./5_pub.outs/2_supp_tables/Figure_S6_1_DEG_enrich.xlsx"
plotDF <- read.xlsx(fn)


col1 <- c("Bcell"="#4daf4a", "Monocyte"="#984ea3", "NKcell"="#aa4b56", "Tcell"="#ffaa00")
col2 <- c("LPS"="#fb9a99", "LPS+DEX"="#e31a1c", "PHA"="#a6cee3", "PHA+DEX"="#1f78b4")


plotDF <- plotDF%>%
    mutate(condition2=gsub("-", "+", conditions), MCls=gsub("_.*", "", conditions),
           contrast=gsub("-", "+", gsub(".*_", "", conditions)),
           log_odds=log(odds), log_lower=log(CI_lower), log_upper=log(CI_upper))

p <- ggplot(plotDF, aes(x=log_odds, y=condition2))+
   geom_errorbarh(aes(xmax=log_upper, xmin=log_lower, colour=contrast),
       size=0.5, height=0.2)+ 
   geom_point(aes(colour=contrast), shape=19, size=1.5)+
   scale_colour_manual(values=col2)+
   geom_vline(aes(xintercept=0), size=0.25, linetype="dashed")+ 
   xlab("log odds ratio")+xlim(-3, 1)+
   ###scale_y_discrete(labels=ylab2)+ 
   theme_bw()+
   theme(##plot.title=element_text(hjust=0.5, size=14),
         axis.title.y=element_blank(),
         axis.title.x=element_text(size=10),
         axis.text.x=element_text(size=10),
         axis.text.y=element_text(size=10),
         ## legend.position="none")
         legend.title=element_blank(),
         legend.text=element_text(size=9),
         legend.key.size=unit(0.4, "cm"))
         ## legend.position="none")


figfn <- paste(outdir, "FigS5_3_olapDEG_enrich.forest.pdf", sep="")
pdf(figfn, width=5.5, height=5)
print(p)
dev.off()

## figfn <- paste(outdir, "Figure_S6_1_DEG_enrich.forest.png", sep="")
## ggsave(figfn, p, device="png", width=550, height=480, units="px", dpi=120)    





#############
##############





############################
### enrichment analysis ####
#############################



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

