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
library(openxlsx)
library(cowplot)
library(ggrepel)
##

rm(list=ls())

###
###
outdir <- "./3_twas_summary.outs/"
if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)


###########################
### calculate FDR 
###########################

trait <- sort(read.table("traits.txt")$V1)[3]

methods <- c("aloft_minP", "aloft_topPIP_union", "gtex_minP", "gtex_topPIP_union")


for ( ii in methods){
   ##
   cat(ii, "\n") 
   fn <- paste("./twas_smr.outs/", trait, "_", ii, "_twas.txt.gz", sep="")
   res <- fread(fn, header=T, data.table=F) 
   res <- res%>%mutate(FDR=p.adjust(pval_gwas, "BH"))
   ## 
   opfn <- paste(outdir, trait, "_", ii, "_twas.txt.gz", sep="")
   fwrite(res, file=opfn, sep="\t", quote=F, na=NA)
   ##
}    

######################
### autosome genes ###
######################

## autosome <- as.character(1:22) 
## grch38_unq <- grch38%>%
##     dplyr::filter(chr%in%autosome, grepl("protein", biotype))%>%
##     distinct(ensgene, chr, .keep_all=T)%>%dplyr::select(gene=ensgene, chr, biotype, symbol)


#########################################
### histogram distribution of p-value ###
#########################################


traits <- sort(read.table("traits.txt")$V1)
trait_short <- c("Asthma_GBMI_afr", "Asthma_GBMI_eur", "Asthma_GBMI_full", "Asthma_UKB_eur")
names(trait_short) <- traits


methods <- c("aloft_minP", "aloft_topPIP_union",
             "gtex_minP", "gtex_topPIP_union")

 
trait <- traits[3]
trait0 <- trait_short[trait]

###    
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
  
p2 <- ggplot()+
   geom_histogram(data=plotDF, aes(x=pval_gwas, color=conditions), fill=NA, bins=40)+
   geom_text(data=annoDF, aes(x=xpos, y=ypos, label=eq), size=3, parse=T)+
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
         axis.text=element_text(size=10),
         axis.title=element_text(size=10),
         strip.text=element_text(size=10))
 
figfn <- paste(outdir, "Figure1_", trait, "_pval.hist.png", sep="")
ggsave(figfn, p2, width=580, height=550, units="px", dpi=120)



############################
### Density plots 
############################



traits <- sort(read.table("traits.txt")$V1)
trait_short <- c("Asthma_GBMI_afr", "Asthma_GBMI_eur", "Asthma_GBMI_full", "Asthma_UKB_eur")
names(trait_short) <- traits


methods <- c("aloft_minP", "aloft_topPIP_union", "gtex_minP", "gtex_topPIP_union")

 
trait <- traits[3]
trait0 <- trait_short[trait]

###    
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



###
###
plotDF2 <- plotDF%>%mutate(gr2=ifelse(grepl("aloft", conditions), "ALOFT", "GTEx"))

 
p0 <- ggplot(plotDF2, aes(x=pval_gwas, colour=conditions))+
   geom_density()+
   scale_color_manual(values=c("aloft_minP"="#f1b6da", "aloft_topPIP_union"="#d01c8b",
                               "gtex_minP"="#b8e186", "gtex_topPIP_union"="#4dac26"),
       labels=c("aloft_minP"="ALOFT_SMR_standard", "aloft_topPIP_union"="ALOFT_SMR_annotation",
                "gtex_minP"="GTEx_SMR_standard", "gtex_topPIP_union"="GTEx_SMR_annotation"))+
   facet_wrap(~gr2, nrow=1, ncol=2, scales="fixed")+ 
   xlab(bquote(~italic(plain(P))~"of TWAS"))+
   ylab("Density #genes")+
   theme_bw()+
   theme(legend.title=element_blank(),
         legend.text=element_text(size=8),
         legend.key.size=grid::unit(0.6, "lines"),
         legend.box.background=element_blank(),
         legend.background=element_blank(),
         axis.text=element_text(size=10),
         axis.title=element_text(size=10),
         strip.text=element_text(size=12))

###
###

figfn <- paste(outdir, "Figure1.2_", trait, "_pval.density.png", sep="")
ggsave(figfn, p0, width=850, height=420, units="px", dpi=120)


### aloft test
pval0 <- plotDF2%>%filter(conditions=="aloft_minP")%>%pull(pval_gwas)
pval1 <- plotDF2%>%filter(conditions=="aloft_topPIP_union")%>%pull(pval_gwas)
wilcox.test(pval0, pval1, alternative="greater")

### GTEx test
pval0 <- plotDF2%>%filter(conditions=="gtex_minP")%>%pull(pval_gwas)
pval1 <- plotDF2%>%filter(conditions=="gtex_topPIP_union")%>%pull(pval_gwas)
wilcox.test(pval0, pval1, alternative="greater")




###
### END









######################
### qq plots 
######################

traits <- sort(read.table("traits.txt")$V1)

trait_short <- c("Asthma_GBMI_afr", "Asthma_GBMI_eur", "Asthma_GBMI_full", "Asthma_UKB_eur")
names(trait_short) <- traits

methods <- c("aloft_minP", "aloft_topPIP_union", "aloft_topPIP_dtss",
             "gtex_minP", "gtex_topPIP_union", "gtex_topPIP_dtss")
col2 <- c("aloft_minP"="#00BFC4", "aloft_topPIP_union"="#F8766D", "aloft_topPIP_dtss"="#7CAE00",
             "gtex_minP"="#00BFC4", "gtex_topPIP_union"="#F8766D", "gtex_topPIP_dtss"="#7CAE00")

for ( trait in traits[3:4]){

 
   plotDF <- map_dfr(methods[1:3], function(ii){
      ##
      fn <- paste("./twas_smr.outs/", trait, "_", ii, "_twas.txt.gz", sep="")
      res0 <- fread(fn, header=T,  data.table=F)
      res0 <- res0[,c("gene", "pval_gwas")]
      ngene <- nrow(res0) 
      res0 <- res0%>%arrange(pval_gwas)%>%
          mutate(observed=-log10(pval_gwas), expected=-log10(ppoints(ngene)), conditions=ii)
      ##
      res0
   })

   trait0 <- trait_short[trait]
    
   p1  <- ggplot(plotDF, aes(x=expected, y=observed, color=conditions))+
      geom_point(size=0.1)+
      scale_color_manual(values=col2, guide=guide_legend(override.aes=list(size=1.5)))+ 
   ## geom_ribbon(aes(ymax=cupper, ymin=clower), fill="grey", alpha=0.5)+ 
      geom_abline(slope=1, intercept=0)+
      xlab(bquote("Expected"~ -log[10]~"("~italic(plain(P))~")"))+
      ylab(bquote("Observed"~-log[10]~"("~italic(plain(P))~")"))+
      ggtitle(paste(trait0, "(ALOFT)"))+ 
      theme_bw()+
      theme(
         legend.position=c(0.3, 0.8),  
         legend.title=element_blank(),
         legend.background=element_blank(),
         legend.box.background=element_blank(),
         legend.key.size=grid::unit(0.6, "lines"),
         legend.text=element_text(size=8),
         axis.text=element_text(size=10),
         axis.title=element_text(size=10),
         plot.title=element_text(hjust=0.5, size=10))
  ## ,
  ##        panel.grid.major=element_blank(),
  ##        panel.grid.minor=element_blank())
    
   ### GTEX
   plotDF2 <- map_dfr(methods[4:6], function(ii){
      ##
      fn <- paste("./twas_smr.outs/", trait, "_", ii, "_twas.txt.gz", sep="")
      res0 <- fread(fn, header=T,  data.table=F)
      res0 <- res0[,c("gene", "pval_gwas")]
      ngene <- nrow(res0) 
      res0 <- res0%>%arrange(pval_gwas)%>%
          mutate(observed=-log10(pval_gwas), expected=-log10(ppoints(ngene)), conditions=ii)
      ##
      res0
   })
   
   p2  <- ggplot(plotDF2, aes(x=expected, y=observed, color=conditions))+
      geom_point(size=0.1)+
      scale_color_manual(values=col2, guide=guide_legend(override.aes=list(size=1.5)))+ 
   ## geom_ribbon(aes(ymax=cupper, ymin=clower), fill="grey", alpha=0.5)+ 
      geom_abline(slope=1, intercept=0)+
      xlab(bquote("Expected"~ -log[10]~"("~italic(plain(P))~")"))+
      ylab(bquote("Observed"~-log[10]~"("~italic(plain(P))~")"))+
      ggtitle(paste(trait0, "(GTEx)"))+ 
      theme_bw()+
      theme(
         legend.position=c(0.3, 0.8), 
         legend.title=element_blank(),
         legend.background=element_blank(),
         legend.box.background=element_blank(),
         legend.key.size=grid::unit(0.6, "lines"),
         legend.text=element_text(size=8),
         axis.text=element_text(size=10),
         axis.title=element_text(size=10),
         plot.title=element_text(hjust=0.5, size=10))
        ## ,
        ##  panel.grid.major=element_blank(),
        ##  panel.grid.minor=element_blank())


     pcomb <- plot_grid(p1, p2, nrow=1, ncol=2, align="h", axis="tb")
     figfn <- paste(outdir, "Figure2_", trait, ".qq.png", sep="")
     ggsave(figfn, pcomb, width=780, height=380, units="px", dpi=120)
  
     ###
     cat(trait0, "\n")    
}    
       




##################################################
### scatter plots, topPIP_union vs topPIP_dtss 
###################################################


## outdir2 <- "./3_twas_summary.outs/3.0_compare_annot/"
## if ( !file.exists(outdir2)) dir.create(outdir2, showWarnings=F, recursive=T)


## for (trait in traits[3:4]){

##     ##
##     trait0 <- trait_short[trait]
    
##     ################# 
##     #### ALOFT
##     ##################
##     ### aloft_topPIP_union
##     ii <- methods[2]
##     fn <- paste("./twas_smr.outs/", trait, "_", ii, "_twas.txt.gz", sep="")
##     res0 <- fread(fn, header=T,  data.table=F)
##     res0 <- res0[,c("gene", "pval_gwas")] 
##     res0 <- res0%>%
##         mutate(qval=qvalue(pval_gwas)$qvalues, is_sig_y=ifelse(qval<0.1, 1, 0),
##                log10p_y=-log10(pval_gwas))%>%
##         dplyr::select(gene, log10p_y, is_sig_y)
    
    
    
##     ### aloft_topPIP_dtss
##     ii <- methods[3]
##     fn <- paste("./twas_smr.outs/", trait, "_", ii, "_twas.txt.gz", sep="")
##     res2 <- fread(fn, header=T,  data.table=F)
##     res2 <- res2[,c("gene", "pval_gwas")]
##     res2 <- res2%>%
##         mutate(qval=qvalue(pval_gwas)$qvalues, is_sig_x=ifelse(qval<0.1, 1, 0),
##                log10p_x=-log10(pval_gwas))%>%
##         dplyr::select(gene, log10p_x, is_sig_x)

##     plotDF <- res0%>%inner_join(res2, by="gene")%>%
##         mutate(sig_gr=case_when(is_sig_y==1&is_sig_x==0~"gr1",
##                                 is_sig_y==0&is_sig_x==1~"gr2",
##                                 is_sig_y==1&is_sig_x==1~"gr3",
##                                 TRUE~"gr0"))
    
##     ###
##     p1 <- ggplot(plotDF, aes(x=log10p_x, y=log10p_y, color=sig_gr))+
##         geom_point(size=0.8)+
##        geom_abline(slope=1, intercept=0, color="blue", linewidth=0.4)+
##        scale_color_manual(values=c("gr0"="grey", "gr1"="red", "gr2"="blue", "gr3"="green"),
##            labels=c("gr1"="sig in peak", "gr2"="sig in dtss", "gr3"="both"),
##            limits=c("gr1", "gr2", "gr3"), guide=guide_legend(override.aes=list(size=1.2)))+         
##        xlab(bquote("DTSS"~-log[10]~"("~italic(plain(P))~")"))+
##       ylab(bquote("Peak"~-log[10]~"("~italic(plain(P))~")"))+
##       ggtitle(paste(trait0, "(ALOFT)"))+ 
##       theme_bw()+
##       theme(legend.position="none",
##          legend.title=element_blank(),
##          legend.background=element_blank(),
##          legend.box.background=element_blank(),
##          legend.key.size=grid::unit(0.6, "lines"),
##          legend.text=element_text(size=8),
##           axis.text=element_text(size=10),
##          axis.title=element_text(size=10),
##          plot.title=element_text(hjust=0.5, size=10))
    

##     ##############
##     ### GTEx 
##     ###############
    
##     ###
##     ### GTEx_topPIP_union
##     ii <- methods[5]
##     fn <- paste("./twas_smr.outs/", trait, "_", ii, "_twas.txt.gz", sep="")
##     res0 <- fread(fn, header=T,  data.table=F)
##     res0 <- res0[,c("gene", "pval_gwas")] 
##     res0 <- res0%>%
##         mutate(qval=qvalue(pval_gwas)$qvalues, is_sig_y=ifelse(qval<0.1, 1, 0),
##                log10p_y=-log10(pval_gwas))%>%
##         dplyr::select(gene, log10p_y, is_sig_y)
    
    
    
##     ### GTEx_topPIP_dtss
##     ii <- methods[6]
##     fn <- paste("./twas_smr.outs/", trait, "_", ii, "_twas.txt.gz", sep="")
##     res2 <- fread(fn, header=T,  data.table=F)
##     res2 <- res2[,c("gene", "pval_gwas")]
##     res2 <- res2%>%
##         mutate(qval=qvalue(pval_gwas)$qvalues, is_sig_x=ifelse(qval<0.1, 1, 0),
##                log10p_x=-log10(pval_gwas))%>%
##         dplyr::select(gene, log10p_x, is_sig_x)

##     plotDF2 <- res0%>%inner_join(res2, by="gene")%>%
##         mutate(sig_gr=case_when(is_sig_y==1&is_sig_x==0~"gr1",
##                                 is_sig_y==0&is_sig_x==1~"gr2",
##                                 is_sig_y==1&is_sig_x==1~"gr3",
##                                 TRUE~"gr0"))
    
##     ###
##     p2 <- ggplot(plotDF2, aes(x=log10p_x, y=log10p_y, color=sig_gr))+
##         geom_point(size=0.8)+
##        geom_abline(slope=1, intercept=0, color="blue", linewidth=0.4)+
##        scale_color_manual(values=c("gr0"="grey", "gr1"="red", "gr2"="blue", "gr3"="green"),
##            labels=c("gr1"="sig in peak", "gr2"="sig in dtss", "gr3"="both"),
##            limits=c("gr1", "gr2", "gr3"), guide=guide_legend(override.aes=list(size=1.2)))+         
##        xlab(bquote("DTSS"~-log[10]~"("~italic(plain(P))~")"))+
##       ylab(bquote("Peak"~-log[10]~"("~italic(plain(P))~")"))+
##       ggtitle(paste(trait0, "(GTEx)"))+ 
##       theme_bw()+
##       theme(legend.title=element_blank(),
##          legend.background=element_blank(),
##          legend.box.background=element_blank(),
##          legend.key.size=grid::unit(0.6, "lines"),
##          legend.text=element_text(size=8),
##           axis.text=element_text(size=10),
##          axis.title=element_text(size=10),
##          plot.title=element_text(hjust=0.5, size=10))    

   
##     ###
##     ###
##      pcomb <- plot_grid(p1, p2, nrow=1, ncol=2, rel_widths=c(1, 1.4), align="h", axis="tb")
##      figfn <- paste(outdir2, "Figure1.1_", trait, "_compare_annot_log10p_.scatter.png", sep="")
##      ggsave(figfn, pcomb, width=780, height=380, units="px", dpi=120)

## }



###########################
### compare zscore
#########################
 
## for (trait in traits[3:4]){

##     ##
##     trait0 <- trait_short[trait]
    
##     ################# 
##     #### ALOFT
##     ##################
##     ### aloft_topPIP_union
##     ii <- methods[2]
##     fn <- paste("./twas_smr.outs/", trait, "_", ii, "_twas.txt.gz", sep="")
##     res0 <- fread(fn, header=T,  data.table=F)
##     res0 <- res0[,c("gene", "pval_gwas", "zscore_gwas")] 
##     res0 <- res0%>%
##         mutate(qval=qvalue(pval_gwas)$qvalues, is_sig_y=ifelse(qval<0.1, 1, 0),
##                log10p_y=-log10(pval_gwas))%>%
##         dplyr::select(gene, log10p_y, is_sig_y, zscore_y=zscore_gwas)
    
    
    
##     ### aloft_topPIP_dtss
##     ii <- methods[3]
##     fn <- paste("./twas_smr.outs/", trait, "_", ii, "_twas.txt.gz", sep="")
##     res2 <- fread(fn, header=T,  data.table=F)
##     res2 <- res2[,c("gene", "pval_gwas", "zscore_gwas")]
##     res2 <- res2%>%
##         mutate(qval=qvalue(pval_gwas)$qvalues, is_sig_x=ifelse(qval<0.1, 1, 0),
##                log10p_x=-log10(pval_gwas))%>%
##         dplyr::select(gene, log10p_x, is_sig_x, zscore_x=zscore_gwas)

##     plotDF <- res0%>%inner_join(res2, by="gene")%>%
##         mutate(sig_gr=case_when(is_sig_y==1&is_sig_x==0~"gr1",
##                                 is_sig_y==0&is_sig_x==1~"gr2",
##                                 is_sig_y==1&is_sig_x==1~"gr3",
##                                 TRUE~"gr0"))
    
##     ###
##     p1 <- ggplot(plotDF, aes(x=zscore_x, y=zscore_y, color=sig_gr))+
##         geom_point(size=0.8)+
##        geom_abline(slope=1, intercept=0, color="blue", linewidth=0.4)+
##        scale_color_manual(values=c("gr0"="grey", "gr1"="red", "gr2"="blue", "gr3"="green"),
##            labels=c("gr1"="sig in peak", "gr2"="sig in dtss", "gr3"="both"),
##            limits=c("gr1", "gr2", "gr3"), guide=guide_legend(override.aes=list(size=1.2)))+         
##       xlab(bquote(italic(Z)~"-score from DTSS"))+
##       ylab(bquote(italic(Z)~"-score from Peak"))+
##       ggtitle(paste(trait0, "(ALOFT)"))+ 
##       theme_bw()+
##       theme(legend.position="none",
##          legend.title=element_blank(),
##          legend.background=element_blank(),
##          legend.box.background=element_blank(),
##          legend.key.size=grid::unit(0.6, "lines"),
##          legend.text=element_text(size=8),
##           axis.text=element_text(size=10),
##          axis.title=element_text(size=10),
##          plot.title=element_text(hjust=0.5, size=10))
    

##     ##############
##     ### GTEx 
##     ###############
    
##     ###
##     ### GTEx_topPIP_union
##     ii <- methods[5]
##     fn <- paste("./twas_smr.outs/", trait, "_", ii, "_twas.txt.gz", sep="")
##     res0 <- fread(fn, header=T,  data.table=F)
##     res0 <- res0[,c("gene", "pval_gwas", "zscore_gwas")] 
##     res0 <- res0%>%
##         mutate(qval=qvalue(pval_gwas)$qvalues, is_sig_y=ifelse(qval<0.1, 1, 0),
##                log10p_y=-log10(pval_gwas))%>%
##         dplyr::select(gene, zscore_y=zscore_gwas, is_sig_y)
    
    
    
##     ### GTEx_topPIP_dtss
##     ii <- methods[6]
##     fn <- paste("./twas_smr.outs/", trait, "_", ii, "_twas.txt.gz", sep="")
##     res2 <- fread(fn, header=T,  data.table=F)
##     res2 <- res2[,c("gene", "pval_gwas", "zscore_gwas")]
##     res2 <- res2%>%
##         mutate(qval=qvalue(pval_gwas)$qvalues, is_sig_x=ifelse(qval<0.1, 1, 0),
##                log10p_x=-log10(pval_gwas))%>%
##         dplyr::select(gene, zscore_x=zscore_gwas, is_sig_x)

##     plotDF2 <- res0%>%inner_join(res2, by="gene")%>%
##         mutate(sig_gr=case_when(is_sig_y==1&is_sig_x==0~"gr1",
##                                 is_sig_y==0&is_sig_x==1~"gr2",
##                                 is_sig_y==1&is_sig_x==1~"gr3",
##                                 TRUE~"gr0"))
    
##     ###
##     p2 <- ggplot(plotDF2, aes(x=zscore_x, y=zscore_y, color=sig_gr))+
##         geom_point(size=0.8)+
##        geom_abline(slope=1, intercept=0, color="blue", linewidth=0.4)+
##        scale_color_manual(values=c("gr0"="grey", "gr1"="red", "gr2"="blue", "gr3"="green"),
##            labels=c("gr1"="sig in peak", "gr2"="sig in dtss", "gr3"="both"),
##            limits=c("gr1", "gr2", "gr3"), guide=guide_legend(override.aes=list(size=1.2)))+
##       xlab(bquote(italic(Z)~"-score from DTSS"))+
##       ylab(bquote(italic(Z)~"-score from Peak"))+        
##       ggtitle(paste(trait0, "(GTEx)"))+ 
##       theme_bw()+
##       theme(legend.title=element_blank(),
##          legend.background=element_blank(),
##          legend.box.background=element_blank(),
##          legend.key.size=grid::unit(0.6, "lines"),
##          legend.text=element_text(size=8),
##           axis.text=element_text(size=10),
##          axis.title=element_text(size=10),
##          plot.title=element_text(hjust=0.5, size=10))    

   
##     ###
##     ###
##      pcomb <- plot_grid(p1, p2, nrow=1, ncol=2, rel_widths=c(1, 1.4), align="h", axis="tb")
##      figfn <- paste(outdir2, "Figure1.2_", trait, "_compare_annot_zscore_.scatter.png", sep="")
##      ggsave(figfn, pcomb, width=780, height=380, units="px", dpi=120)

## }






#########################################
### union vs fast minimum P
####################################

## for (trait in traits[3:4]){

##     ##
##     trait0 <- trait_short[trait]

##     ###      
##     df0 <- map_dfr(methods[1:2], function(ii){
##        ## 
##        fn <- paste("./twas_smr.outs/", trait, "_", ii, "_twas.txt.gz", sep="")
##        res0 <- fread(fn, header=T,  data.table=F)
##        res0 <- res0[,c("gene", "pval_gwas")]
##        res0 <- res0%>%mutate(log10_pval=-log10(pval_gwas), conditions=ii)%>%dplyr::select(-pval_gwas)
##        res0 
##     })    
##     plotDF <- df0%>%pivot_wider(id_cols=gene, names_from=conditions, values_from=log10_pval)
##     names(plotDF) <- c("gene", "minP", "topPIP")
##     plotDF <- plotDF%>%drop_na(minP, topPIP)
##     ###
##     p1 <- ggplot(plotDF, aes(x=minP, y=topPIP))+
##         geom_point(color="grey", size=0.8)+
##        geom_abline(slope=1, intercept=0, color="blue", linewidth=0.4)+
##        xlab(bquote("fastqtl"~-log[10]~"("~italic(plain(P))~")"))+
##       ylab(bquote("dap-g-peak"~-log[10]~"("~italic(plain(P))~")"))+
##       ggtitle(paste(trait0, "(ALOFT)"))+ 
##       theme_bw()+
##       theme(axis.text=element_text(size=10),
##          axis.title=element_text(size=10),
##          plot.title=element_text(hjust=0.5, size=10))    


##      ###
##      ###      
##     df0 <- map_dfr(methods[4:5], function(ii){
##        ## 
##        fn <- paste("./twas_smr.outs/", trait, "_", ii, "_twas.txt.gz", sep="")
##        res0 <- fread(fn, header=T,  data.table=F)
##        res0 <- res0[,c("gene", "pval_gwas")]
##        res0 <- res0%>%mutate(log10_pval=-log10(pval_gwas), conditions=ii)%>%dplyr::select(-pval_gwas)
##        res0 
##     })    
##     plotDF2 <- df0%>%pivot_wider(id_cols="gene", names_from=conditions, values_from=log10_pval)
##     names(plotDF2) <- c("gene", "minP", "topPIP")
##     plotDF2 <- plotDF2%>%drop_na(minP, topPIP)
     
##     p2 <- ggplot(plotDF2, aes(x=minP, y=topPIP))+
##         geom_point(color="grey", size=0.8)+
##        geom_abline(slope=1, intercept=0, color="blue", linewidth=0.4)+
##        xlab(bquote("fastqtl"~-log[10]~"("~italic(plain(P))~")"))+
##       ylab(bquote("dap-g-peak"~-log[10]~"("~italic(plain(P))~")"))+
##       ggtitle(paste(trait0, "(GTEx)"))+ 
##       theme_bw()+
##       theme(axis.text=element_text(size=10),
##          axis.title=element_text(size=10),
##          plot.title=element_text(hjust=0.5, size=10))     
    
##     ###
##     ###
##      pcomb <- plot_grid(p1, p2, nrow=1, ncol=2, align="h", axis="tb")
##      figfn <- paste(outdir, "Figure3.2_", trait, "_fastqtl.scatter.png", sep="")
##      ggsave(figfn, pcomb, width=780, height=380, units="px", dpi=120)

## }





