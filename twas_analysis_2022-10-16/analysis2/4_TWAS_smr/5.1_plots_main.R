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





############################
### Density plots 
############################



trait <- sort(read.table("traits.txt")$V1)[3]
methods <- c("aloft_minP", "aloft_topPIP_union", "gtex_minP", "gtex_topPIP_union")

 
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
### significance test

### aloft test
pval0 <- plotDF2%>%filter(conditions=="aloft_minP")%>%pull(pval_gwas)
pval1 <- plotDF2%>%filter(conditions=="aloft_topPIP_union")%>%pull(pval_gwas)
wil_test <- wilcox.test(pval0, pval1, alternative="greater")
symb <- "Wilcoxon test ***"

### GTEx test
pval0 <- plotDF2%>%filter(conditions=="gtex_minP")%>%pull(pval_gwas)
pval1 <- plotDF2%>%filter(conditions=="gtex_topPIP_union")%>%pull(pval_gwas)
wil_test2 <- wilcox.test(pval0, pval1, alternative="greater")
symb2 <- "Wilcoxon test ***"

anno_df <- data.frame(xpos=c(0.25, 0.25), ypos=c(0.5, 0.5), gr2=c("ALOFT", "GTEx"), symbols=c(symb, symb2))


###
###
plotDF2 <- plotDF%>%mutate(gr2=ifelse(grepl("aloft", conditions), "ALOFT", "GTEx"))

facet_lab2 <- as_labeller(c("ALOFT"="ALOFT bulk", "GTEx"="GTEx WBL"))
 
p0 <- ggplot(plotDF2, aes(x=pval_gwas))+
   geom_density(aes(colour=conditions))+
   scale_color_manual(values=c("aloft_minP"="#f1b6da", "aloft_topPIP_union"="#d01c8b",
                               "gtex_minP"="#b8e186", "gtex_topPIP_union"="#4dac26"),
       labels=c("aloft_minP"="ALOFT_SMR_standard", "aloft_topPIP_union"="ALOFT_SMR_annotation",
                "gtex_minP"="GTEx_SMR_standard", "gtex_topPIP_union"="GTEx_SMR_annotation"))+
   facet_wrap(~gr2, nrow=1, ncol=2, scales="fixed", labeller=facet_lab2)+
   geom_text(data=anno_df, aes(x=xpos, y=ypos, label=symbols), size=3.5)+ 
   xlab("p value of TWAS")+
   ylab("Density #genes")+
   theme_bw()+
   theme(legend.title=element_blank(),
         legend.text=element_text(size=8),
         legend.key.size=grid::unit(0.6, "lines"),
         legend.box.background=element_blank(),
         legend.background=element_blank(),
         legend.position=c(0.8, 0.8),
         axis.text=element_text(size=10),
         axis.title=element_text(size=10),
         strip.text=element_text(size=12))

###
###

figfn <- paste(outdir, "Figure1.2_", trait, "_pval.density.png", sep="")
ggsave(figfn, p0, width=650, height=420, units="px", dpi=120)



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
fn <- paste("./4_INTACT.outs/", trait, "_aloft_topPIP_union_combinfor.txt", sep="")
res <- read.table(fn, header=T, sep="\t")
res <- res%>%left_join(grch38_unq, by="Gene")

res <- res%>%mutate(log10p=-log10(pval_gwas), gr2=factor(chr%%2))

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
   summarize(midpoint=(max(BPcum)+min(BPcum))/2)



sig <- res%>%filter(FDR_twas<0.1)%>%pull(log10p)%>%min()


x <- res%>%dplyr::filter(FDR_twas<0.05)%>%arrange(pval_gwas)
x1 <- x[1:20,]
x2 <- x%>%dplyr::filter(FDR<0.05, peaking_d==3)
annoDF2 <- rbind(x1,x2)%>%distinct(symbol, .keep_all=T)


ymax <- max(res$log10p)+8
p1 <- ggplot(res)+
   geom_point(aes(x=BPcum, y=log10p, color=factor(gr2), size=log10p))+
   geom_text_repel(data=annoDF2,
       aes(x=BPcum, y=log10p, label=symbol), box.padding=0.2,
            max.overlaps=15, fontface="italic", size=3)+      
   geom_hline(yintercept=sig, color="red", linetype="dashed")+ 
   scale_x_continuous("chromosome", label=axis.set$chr, breaks=axis.set$midpoint,
      limits=c(min(res$BPcum), max(res$BPcum)) )+
   scale_y_continuous(bquote(-log[10]~"(p value)"), limits=c(0,ymax))+
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
plotDF <- res%>%dplyr::select(Gene, pval=pval_gwas)%>%arrange(pval)%>%
    mutate(observed=-log10(pval), expected=-log10(ppoints(ngene)),
           clower=-log10(qbeta(p=(1-ci)/2, shape1=seq(ngene), shape2=rev(seq(ngene)))),
           cupper=-log10(qbeta(p=(1+ci)/2, shape1=seq(ngene), shape2=rev(seq(ngene)))) )
           
p2  <- ggplot(plotDF, aes(x=expected, y=observed))+
   geom_point(color="#d01c8b", size=1)+
   ## geom_ribbon(aes(ymax=cupper, ymin=clower), fill="grey", alpha=0.5)+ 
   geom_abline(slope=1, intercept=0)+
   xlab(bquote("Expected"~ -log[10]~"(p value)"))+
   ylab(bquote("Observed"~-log[10]~"(p value)"))+
   theme_bw()+
   theme(legend.position="none",
         axis.text=element_text(size=12),
         axis.title=element_text(size=12),
         panel.grid.major=element_blank(),
         panel.grid.minor=element_blank())

## p2 <- as.grob(p2)


plot_comb <- plot_grid(p1, p2, rel_widths=c(1.7, 1), nrow=1, ncol=2, aligh="h", axis="tb")


###
figfn <- paste(outdir, "Fig1.4_asthma_manhattan_qq_comb.png", sep="")
png(figfn, width=1200, height=450, pointsize=12, res=120)
print(plot_comb)
dev.off()


####################################
### two cohorts Manhanttan plots ### 
####################################

## rm(list=ls())


## autosome <- as.character(1:22) 
## grch38_unq <- grch38%>%
##     dplyr::filter(chr%in%autosome, grepl("protein", biotype))%>%
##     distinct(ensgene, chr, .keep_all=T)%>%dplyr::select(Gene=ensgene, start)



## ### ALOFT 
## fn <- "./5_pub.outs/2_supp_tables/TableS5_1_asthma-risk-genes_ALOFT.txt.gz"
## res <- fread(fn, header=T, data.table=F, sep="\t")
## res <- res%>%dplyr::select(Gene, symbol, chr=chromosome, PIP, peaking_d, pval_twas, FDR_twas, FDR_intact)
 
## res <- res%>%left_join(grch38_unq, by="Gene")
## res <- res%>%mutate(log10p=-log10(pval_twas), gr2=factor(chr%%2))
## res$dataset <- "aloft"


## ### 
## ### add x-axis information 
## ngene <- nrow(res)
## nchr <- max(res$chr)

## ###
## res$BPcum <- NA
## s <- 0
## for (i in 1:nchr){
##    x <- res[res$chr==i, "start"]+s
##    res[res$chr==i,"BPcum"] <- x
##    s <- max(x)
## }




## ############
## ### GTEx ###
## ############
 
## fn <- "./5_pub.outs/2_supp_tables/TableS6_1_asthma-risk-genes_gtex.txt.gz"
## res_gtex <- fread(fn, header=T, data.table=F, sep="\t")
## res_gtex <- res_gtex%>%dplyr::select(Gene, symbol, chr=chromosome, PIP, peaking_d, pval_twas, FDR_twas, FDR_intact)


## res_gtex <- res_gtex%>%left_join(grch38_unq, by="Gene")
## res_gtex <- res_gtex%>%mutate(log10p=log10(pval_twas), gr2=factor(chr%%2))
## res_gtex$dataset <- "gtex"


## ### 
## ### add x-axis information 
## ngene <- nrow(res_gtex)
## nchr <- max(res_gtex$chr)

## ###
## res_gtex$BPcum <- NA
## s <- 0
## for (i in 1:nchr){
##    x <- res_gtex[res_gtex$chr==i, "start"]+s
##    res_gtex[res_gtex$chr==i,"BPcum"] <- x
##    s <- max(x)
## }

## ### plot data
## res_comb <- rbind(res, res_gtex)


## axis.set <- res_comb%>%
##    group_by(chr)%>%
##    summarize(midpoint=(max(BPcum)+min(BPcum))/2)%>%ungroup()



## sig1 <- res_comb%>%filter(dataset=="aloft", FDR_twas<0.1)%>%pull(log10p)%>%min()
## sig0 <- res_comb%>%filter(dataset=="gtex", FDR_twas<0.1)%>%pull(log10p)%>%max()


## ##################################
## ### tables of annotation genes ###
## ##################################

## x <- res_comb%>%dplyr::filter(dataset=="aloft", FDR_twas<0.05, FDR_intact<0.05, peaking_d==3, PIP>0.1)%>%
##     arrange(pval_twas)


## x2 <- res_comb%>%dplyr::filter(dataset=="gtex", FDR_twas<0.05, FDR_intact<0.05, peaking_d==3, PIP>0.1)%>%
##     arrange(pval_twas)

## annoDF2 <- rbind(x, x2)




## ####
## #### manhattan plots 
## ymax <- max(res_comb$log10p)+10
## ymin <- min(res_comb$log10p)-10

## mybreak <- seq(-100, 100, by=25)
## mylab <- abs(mybreak)


## p1 <- ggplot(res_comb)+
##    geom_point(aes(x=BPcum, y=log10p, color=factor(gr2), size=log10p))+
##    geom_text_repel(data=annoDF2,
##        aes(x=BPcum, y=log10p, label=symbol),
##        box.padding=0.2, ##min.segment.length=0, segment.curvature=1e-20,
##        max.overlaps=30, fontface="italic", size=2.5)+      
##    geom_hline(yintercept=c(sig0, 0, sig1),
##               color=c("red", "grey", "red"), linetype=c("dashed", "solid", "dashed"),
##               linewidth=c(0.8, 0.8, 0.8))+ 
##    scale_x_continuous("chromosome", label=axis.set$chr, breaks=axis.set$midpoint,
##       limits=c(min(res$BPcum), max(res$BPcum)) )+
##    scale_y_continuous(bquote(-log[10]~"("~italic(p)~")"), limits=c(ymin, ymax), breaks=mybreak, labels=mylab)+
##    scale_color_manual(values=c("0"="#2171b5", "1"="#6baed6"))+
##    scale_size_continuous(range=c(0.1, 1.5)) +
##    theme_bw()+
##    theme(legend.position="none",
##          axis.title=element_text(size=12),
##          axis.text.x=element_text(size=10, angle=90, vjust=0.5, hjust=1),
##          axis.text.y=element_text(size=12),
##          panel.grid.major=element_blank(),
##          panel.grid.minor=element_blank())



## ################
## ### qq plots ###
## ################

## plotDF <- map_dfr(c("aloft", "gtex"), function(ii){
##     ###
##     res0 <- res_comb%>%filter(dataset==ii)%>%dplyr::select(Gene, pval=pval_twas)%>%arrange(pval)
##     ngene <- nrow(res0)
##     ##
##     res0 <- res0%>%mutate(observed=-log10(pval), expected=-log10(ppoints(ngene)))
##     res0$dataset <- ii
##     res0
## })

## p2  <- ggplot(plotDF, aes(x=expected, y=observed, color=factor(dataset) ))+
##    geom_point(size=1)+
##    scale_color_manual(values=c("aloft"="#d01c8b", "gtex"="#4dac26"),
##        labels=c("aloft"="ALOFT", "gtex"="GTEx"))+    
##    geom_abline(slope=1, intercept=0)+
##    xlab(bquote("Expected"~ -log[10]~"("~italic(plain(P))~")"))+
##    ylab(bquote("Observed"~-log[10]~"("~italic(plain(P))~")"))+
##    theme_bw()+
##    theme(legend.position=c(0.3, 0.8),
##          legend.title=element_blank(),
##          legend.text=element_text(size=8),
##          legend.key.size=grid::unit(0.8, "lines"),
##          legend.box.background=element_blank(),
##          legend.background=element_blank(),         
##          axis.text=element_text(size=12),
##          axis.title=element_text(size=12),
##          panel.grid.major=element_blank(),
##          panel.grid.minor=element_blank())

## ## p2 <- as.grob(p2)


## plot_comb <- plot_grid(p1, p2, rel_widths=c(1.8, 1), nrow=1, ncol=2, aligh="h", axis="tb")


## ###
## figfn <- paste(outdir, "Fig2_asthma_manhattan_qq_comb.png", sep="")
## png(figfn, width=1200, height=450, pointsize=12, res=120)
## print(plot_comb)
## dev.off()










####################################
### summary number of risk genes ###
####################################



trait <- sort(read.table("traits.txt")$V1)[3]
methods <- c("aloft_minP", "aloft_topPIP_union", "gtex_minP", "gtex_topPIP_union")


###
##
## autosome <- as.character(1:22) 
## grch38_unq <- grch38%>%
##     dplyr::filter(chr%in%autosome, grepl("protein", biotype))%>%
##     distinct(ensgene, chr, .keep_all=T)%>%dplyr::select(gene=ensgene, chr, biotype, symbol)


###
### aloft_min
ii <- methods[1]
fn0 <- paste("./3_twas_summary.outs/", trait, "_", ii, "_twas.txt.gz", sep="")
res0 <- fread(fn0, header=T, data.table=F)
aloft0 <- res0%>%filter(FDR<0.1)%>%pull(gene)%>%unique()


###
### aloft_topPIP
ii <- methods[2]
fn <- paste("./3_twas_summary.outs/", trait, "_", ii, "_twas.txt.gz", sep="")
res <- fread(fn, header=T, data.table=F)
aloft2 <- res%>%filter(FDR<0.1)%>%pull(gene)%>%unique()



###
### gtex_min
ii <- methods[3]
fn0 <- paste("./3_twas_summary.outs/", trait, "_", ii, "_twas.txt.gz", sep="")
res0 <- fread(fn0, header=T, data.table=F)
gtex0<- res0%>%filter(FDR<0.1)%>%pull(gene)%>%unique()


###
### gtex_topPIP
ii <- methods[4]
fn <- paste("./3_twas_summary.outs/", trait, "_", ii, "_twas.txt.gz", sep="")
res <- fread(fn, header=T, data.table=F)
gtex2 <- res%>%filter(FDR<0.1)%>%pull(gene)%>%unique()

(length(gtex2)-length(gtex0))/length(gtex0)


sig0 <- union(aloft0, gtex0)
sig2 <- union(aloft2, gtex2)




## nsig <- sum(res$FDR<0.1)
ptwas <- read.table("PTWAS.gtex_v8.all_traits.tgenes.fdr_05.txt", header=F)
names(ptwas) <- c("traits", "gene", "tissue", "FDR") 
ptwas <- ptwas%>%mutate(gene2=gsub("\\..*", "", gene))
ptwas2 <- ptwas%>%filter(grepl("asthma|Asthma", traits)) ##, gene2%in%grch38_unq$gene)

old <- unique(ptwas2$gene2)

olap <- intersect(gene_asthma, old) ## 359 gene




##############################
### summary INTACT results ###
##############################


trait <- sort(read.table("traits.txt")$V1)[3]
methods <- c("aloft_minP", "aloft_topPIP_union", "gtex_minP", "gtex_topPIP_union")

### minP
## aloft
fn <- paste("./4_INTACT.outs/", trait, "_aloft_minP_intact.txt", sep="")
res0 <- fread(fn, header=T)
aloft0 <- res0%>%filter(FDR<0.1)%>%pull(Gene)%>%unique()

##gtex
fn <- paste("./4_INTACT.outs/", trait, "_gtex_minP_intact.txt", sep="")
res0 <- fread(fn, header=T)
gtex0 <- res0%>%filter(FDR<0.1)%>%pull(Gene)%>%unique()


###topPIP
## aloft
fn <- paste("./4_INTACT.outs/", trait, "_aloft_topPIP_union_intact.txt", sep="")
res2 <- fread(fn, header=T)
aloft2 <- res2%>%filter(FDR<0.1)%>%pull(Gene)%>%unique()

##gtex
fn <- paste("./4_INTACT.outs/", trait, "_gtex_topPIP_union_intact.txt", sep="")
res2 <- fread(fn, header=T)
gtex2 <- res2%>%filter(FDR<0.1)%>%pull(Gene)%>%unique()


sig0 <- union(aloft0, gtex0)
sig2 <- union(aloft2, gtex2)


###
###

fn <- "./5_pub.outs/2_supp_tables/TableS5_1_asthma-risk-genes_ALOFT.txt.gz"
res <- fread(fn, header=T, sep="\t")
sig0 <- res%>%filter(FDR_intact<0.1, peaking_d==3)%>%pull(Gene)%>%unique()


###
fn <- "./5_pub.outs/2_supp_tables/TableS6_1_asthma-risk-genes_gtex.txt.gz"
res2 <- fread(fn, header=T, sep="\t")
sig2 <- res2%>%filter(FDR_intact<0.1, peaking_d==3)%>%pull(Gene)%>%unique()

sig2 <- union(sig0, sig2)
