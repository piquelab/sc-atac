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

library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(viridis)
library(openxlsx)
library(cowplot)
library(ggrepel)
##

options(scipen=16)

rm(list=ls())


outdir <- "./5_pub.outs/2_supp_tables/"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)




#######################
### supp table 
########################

trait <- sort(read.table("traits.txt")$V1)[3]

fn <- paste("./4_INTACT.outs/", trait, "_gtex_topPIP_union_combinfor.txt", sep="")
res <- read.table(fn, header=T, sep="\t")


###
opfn <- gzfile(paste(outdir, "TableS6_1_asthma-risk-genes_gtex.txt.gz", sep=""))
write.table(res, file=opfn, quote=F, row.names=F, sep="\t")






####################################
### annotation genetic variants 
#####################################


fn <- paste(outdir, "TableS6_1_asthma-risk-genes_gtex.txt.gz", sep="")
res <- fread(fn, header=T, data.table=F)

## res <- res[, -5]
## ### output 
## opfn <- gzfile(paste(outdir, "TableS6_1_asthma-risk-genes_gtex.txt.gz", sep=""))
## write.table(res, file=opfn, quote=F, row.names=F, sep="\t")




res2 <- res%>%filter(FDR_intact<0.1, peaking_d==3) 


###
### Add peak infor
fn <- "/nfs/rprdata/julong/sc-atac/analyses.2021-02-05/2_Differential/1.3_DiffPeak.outs/3.0_DESeq_indi.results.rds"
resDAR <- read_rds(fn)%>%mutate(is_sig=ifelse(p.adjusted<0.1&abs(estimate)>0.5, sign(estimate), 0))

##
gene <- unique(resDAR$gene)
cvt <- str_split(gene, "-", simplify=T) 
DF <- data.frame(peak=gene, chr=as.numeric(cvt[,1]), pos1=as.integer(cvt[,2]), pos2=as.integer(cvt[,3]))

##
peakSel <- sapply(1:nrow(res2), function(i){
   ##
   chr_pos <- unlist(str_split(res2$chr_pos_grch37[i], "_", simplify=T)) 
   chri <- as.integer(chr_pos[1])
   posi <- as.integer(chr_pos[2]) 
   ###
   peakSel <- DF%>%filter(chr==chri, pos1<=posi, pos2>=posi)%>%pull(peak)
   peak0 <- paste(peakSel, collapse=";") 
})

res2$peakAnno <- as.character(peakSel)


#############################
### Add response motif
############################

prefix <- "/nfs/rprdata/julong/sc-atac/genetic.analysis_torus_2021-10-14/SNPannotation/"

###
fn <- paste(prefix, "2.3_scanPwm_output/zzz_allmotif.bed.txt.gz", sep="")
snpmotif <- fread(fn, header=F, data.table=F)
names(snpmotif) <- c("motifs", "chr", "pos", "variants", "score_ref", "score_alt")
snpmotif <- snpmotif%>%mutate(chr_pos=paste(chr, pos, sep="_"))


###
### response motifs
contrast <- c("LPS", "LPS-DEX", "PHA", "PHA-DEX")
MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
###
combs <- paste(rep(MCls, each=4), rep(contrast, times=4), sep="_")
##
motifRsps <- map_dfr(combs, function(ii){
   ##
   prefix <- "/nfs/rprdata/julong/sc-atac/analyses.2021-02-05/3_motif/2_motif.activities.outs/MotifList/" 
   fn <- paste(prefix, "1_comb_", ii, "_motif.txt", sep="")
   tmp <- read.table(fn)
   names(tmp) <- "motif.id"
   tmp$conditions <- ii
   tmp
})    

motifRsps <- sort(unique(motifRsps$motif.id))


snpmotif2 <- snpmotif%>%filter(motifs%in%motifRsps)

###
annoMotif <- sapply(1:nrow(res2), function(i){    
  ##   
  chr_pos_i <- res2$chr_pos_grch37[i] 
  ##  
  motifSel <- snpmotif2%>%filter(chr_pos==chr_pos_i)%>%
      mutate(motif_score=paste(motifs, round(score_ref, 3), round(score_alt, 3), sep="_"))%>%
      pull(motif_score)  
  ###   
  motif0 <- paste(motifSel, collapse=";")
  motif0  
})    
 
res2$motifAnnot <- as.character(annoMotif)


#################
### ase 
#################

fn <- paste(outdir, "TableS6_2_asthma-risk-genes_in_response.txt.gz", sep="")
res <- read.table(fn, header=T, sep="\t")
res2 <- res[,1:17]
    
###
fn2 <- "/nfs/rprdata/julong/sc-atac/demux.2021-01-23/asequant/asefinal.txt.gz"
ase <- fread(fn2, header=T, data.table=F)  
ase <- ase%>%mutate(chr_pos=paste(chr, pos, sep="_"))


###
###
ase_summ <- map_dfr(1:nrow(res2), function(i){
   ##
   pos0 <- res2$chr_pos_grch37[i]
   ase0 <- ase%>%filter(chr_pos==pos0)
   if ( nrow(ase0)>0){
    ###
      ase0 <- ase0%>%slice_min(order_by=pval, n=1)
      if ( nrow(ase0)>1){
          ase0 <- ase0%>%slice_max(order_by=abs(beta), n=1, with_ties=F)
      }    
      ase0 <- ase0%>%dplyr::select(ase_comb=comb, ase_beta=beta, ase_pval=pval, ase_fdr=p.adj)
   }else{
      ase0 <- data.frame(ase_comb=NA, ase_beta=NA, ase_pval=NA, ase_fdr=NA)
   }   
   ase0    
})

res2_comb <- cbind(res2, ase_summ)

### txt.gz format
opfn <- gzfile(paste(outdir, "TableS6_2_asthma-risk-genes_in_response.txt.gz", sep=""))
write.table(res2_comb, file=opfn, quote=F, row.names=F, sep="\t")


####
## fn <- paste(outdir, "TableS6_2_asthma-risk-genes_in_response.txt.gz", sep="")
## res <- read.table(fn, header=T, sep="\t")

## res2 <- res[,-5]
## opfn <- gzfile(paste(outdir, "TableS6_2_asthma-risk-genes_in_response.txt.gz", sep=""))
## write.table(res2, file=opfn, quote=F, row.names=F, sep="\t")

## ### save xlsx file 





## opfn <- paste(outdir, "TableS5_2_asthma-risk-genes_in_response.xlsx", sep="")
## write.xlsx(res2, file=opfn)


### aloft
fn <- paste(outdir, "TableS5_1_asthma-risk-genes_ALOFT.txt.gz", sep="")
res <- fread(fn, header=T, data.table=F, sep="\t")
gene <- res%>%filter(FDR_intact<0.1, peaking_d>0)%>%pull(Gene)%>%unique()

## gtex
fn2 <- paste(outdir, "TableS6_1_asthma-risk-genes_gtex.txt.gz", sep="")
res2 <- fread(fn2, header=T, data.table=F, sep="\t")
gene2 <- res2%>%filter(FDR_intact<0.1, peaking_d>0)%>%pull(Gene)%>%unique()
 
gene2 <- union(gene, gene2)
 


##########################
### differential
##########################


fn <- paste(outdir, "TableS6_2_asthma-risk-genes_in_response.txt.gz", sep="")
res <- read.table(fn, header=T, sep="\t")


fn <- "/nfs/rprdata/julong/sc-atac/analyses.2021-02-05/2_Differential/1.3_DiffPeak.outs/3.0_DESeq_indi.results.rds"
resDAR <- read_rds(fn)%>%mutate(is_sig=ifelse(p.adjusted<0.1&abs(estimate)>0.5, sign(estimate), 0))
names(resDAR)[3] <- "peak"

DF <- NULL
DF <- map_dfr(1:nrow(res), function(i){
    ##
    peak0 <- res$peakAnno[i]
    resDAR2 <- resDAR%>%filter(peak==peak0)%>%
        mutate(Gene=res$Gene[i], symbol=res$symbol[i], genetic_variant=res$genetic_variant[i])
    resDAR2 <- resDAR2%>%
        dplyr::select(Gene, symbol, genetic_variant, MCls, contrast, peak, baseMean, estimate, stderror, statistic,
                      p.value, p.adjusted)
    resDAR2
})

###
opfn <- gzfile(paste(outdir, "TableS6_3_asthma-risk-genes_peak.txt.gz", sep=""))
write.table(DF, file=opfn, quote=F, row.names=F, sep="\t")


####
### save xlsx fil
## fn <- paste(outdir, "TableS5_3_asthma-risk-genes_peak.txt.gz", sep="")
## res <- read.table(fn, header=T, sep="\t")

## opfn <- paste(outdir, "TableS5_3_asthma-risk-genes_peak.xlsx", sep="")
## write.xlsx(res, file=opfn)




    
###########################
### differential motif
###########################

fn <- paste(outdir, "TableS6_2_asthma-risk-genes_in_response.txt.gz", sep="")
res <- read.table(fn, header=T, sep="\t")
 
fn <- "/nfs/rprdata/julong/sc-atac/analyses.2021-02-05/3_motif/2_motif.activities.outs/3_motif.diff.results.rds"
resDiff <- read_rds(fn)%>%mutate(conditions=paste(MCls, contrast, sep="_"))

resDiff <- resDiff%>%dplyr::rename(motif_ID=gene, motif_name=motif)


DF <- map_dfr(1:nrow(res), function(i){
    ##
    motif0 <- gsub("_.*", "", unlist(str_split(res$motifAnnot[i], ";")))
    tmp <- resDiff%>%filter(motif_ID%in%motif0)
    tmp <- tmp%>%mutate(Gene=res$Gene[i], symbol=res$symbol[i], genetic_variant=res$genetic_variant[i])
    ###
    tmp <- tmp%>%
        dplyr::select(Gene, symbol, genetic_variant,
        conditions, MCls, contrast, motif_ID, motif_name, beta, stderr, pval, qval)
    tmp
})

opfn <- gzfile(paste(outdir, "TableS6_4_asthma-risk-genes_motif.txt.gz", sep=""))
write.table(DF, file=opfn, quote=F, row.names=F, sep="\t")


###
### save xlsx
## fn <- paste(outdir, "TableS5_4_asthma-risk-genes_motif.txt.gz", sep="")
## res <- read.table(fn, header=T, sep="\t")

## opfn <- paste(outdir, "TableS5_4_asthma-risk-genes_motif.xlsx", sep="")
## write.xlsx(res, file=opfn)







##################################
#### risk genes from GTEx
##################################










#############################################
### summary the results of risk genes  
#############################################

  
fn <- paste(outdir, "TableS5_1_asthma-risk-genes_ALOFT.txt.gz", sep="")
res <- fread(fn, header=T, sep="\t", data.table=F)
res <- res%>%filter(FDR_intact<0.1)

### gtex
fn <- paste(outdir, "TableS6_1_asthma-risk-genes_gtex.txt.gz", sep="")
res2 <- fread(fn, header=T, sep="\t", data.table=F)
res2 <- res2%>%filter(FDR_intact<0.1)



### DEGs
fn <- "/nfs/rprdata/julong/SCAIP/analyses/SCAIP-B1-6_2020.03.23/6_DEG.CelltypeNew_output/Filter2/2_meta.rds"
res_diff <- read_rds(fn)%>%filter(abs(beta)>0.5&qval<0.1)%>%
    mutate(conditions=paste(MCls, contrast, sep="_"))%>%
    dplyr::select(conditions, gene)

 
plotDF <- map_dfr(sort(unique(res_diff$conditions)), function(ii){
    ##
    DEG <- res_diff%>%filter(conditions==ii)%>%pull(gene)
    olap <- intersect(DEG, unique(res2$Gene, res$Gene))
    df0 <- data.frame(conditions=ii, olap=length(olap), ndeg=length(DEG), prop=round(length(olap)/length(DEG), 3))
    df0
})    

plotDF2 <- plotDF%>%filter(olap>0)%>%
   mutate(condition2=gsub("-", "+", conditions))
 
##
p0 <- ggplot(plotDF2, aes(x=condition2, y=olap))+
    geom_bar(stat="identity", fill="#f768a1")+
    coord_flip()+
    ylab("#DEGs")+    
    theme_bw()+
    theme(legend.title=element_blank(),
          axis.title.x=element_text(size=10),
          axis.title.y=element_blank(),
          axis.text=element_text(size=8))

figfn <- paste(outdir, "Figure_S6_1_DEGs.bar.png", sep="")
ggsave(figfn, p0, width=380, height=450, units="px", dpi=120)


###
### proportion of DEGs
plotDF2 <- plotDF%>%filter(olap>0)%>%
   mutate(prop=100*prop, condition2=gsub("-", "+", conditions))
 
##
p1 <- ggplot(plotDF2, aes(x=condition2, y=prop))+
    geom_bar(stat="identity", fill="#f768a1")+
    coord_flip()+
    ylab("Proportion in DEGs(%)")+    
    theme_bw()+
    theme(legend.title=element_blank(),
          axis.title.x=element_text(size=10),
          axis.title.y=element_blank(),
          axis.text=element_text(size=8))

figfn <- paste(outdir, "Figure_S6_1_1_prop_DEGs.bar.png", sep="")
ggsave(figfn, p1, width=380, height=450, units="px", dpi=120)
 

########################
### test enrichment 
##########################


fn <- paste(outdir, "TableS5_1_asthma-risk-genes_ALOFT.txt.gz", sep="")
res <- fread(fn, header=T, sep="\t", data.table=F)
gene <- unique(res$Gene)
res <- res%>%filter(FDR_intact<0.1)

### gtex
fn <- paste(outdir, "TableS6_1_asthma-risk-genes_gtex.txt.gz", sep="")
res2 <- fread(fn, header=T, sep="\t", data.table=F)
gene2 <- unique(res2$Gene)
res2 <- res2%>%filter(FDR_intact<0.1)

## sig <- rbind(res, res2)%>%filter(peaking_d>0)%>%pull(Gene)%>%unique()



gene_test <- union(gene, gene2)
sig <- union(res$Gene, res2$Gene)


### DEGs
fn <- "/nfs/rprdata/julong/SCAIP/analyses/SCAIP-B1-6_2020.03.23/6_DEG.CelltypeNew_output/Filter2/2_meta.rds"
res_diff <- read_rds(fn)%>%  ###filter(gene%in%gene_test)%>%
    mutate(conditions=paste(MCls, contrast, sep="_"),
           is_sig=ifelse(abs(beta)>0.5&qval<0.1, "is_DEG", "not_DEG"),
           is_twas=ifelse(gene%in%sig, "is_asthma", "not_asthma"))


comb <- sort(unique(res_diff$conditions))
plotDF <- NULL
for ( ii in comb){
    ##
    tmp <- res_diff%>%filter(conditions==ii)
    dat2 <- table(tmp$is_sig, tmp$is_twas)
    n0 <- tmp%>%filter(is_sig=="is_DEG")%>%pull(gene)%>%unique()%>%length()
    ##
    if ( dat2[1,1]>0){ 
       fisher <- fisher.test(dat2) ##, alternative="greater") 
       fisher2 <- fisher.test(dat2, alternative="greater")
    
       df0 <- data.frame(conditions=ii, odds=fisher$estimate,
                      CI_lower=fisher$conf.int[1], CI_upper=fisher$conf.int[2],
                      pval=fisher2$p.value, olap=dat2[1,1], ndeg=n0) 
       ##
       plotDF <- rbind(plotDF, df0)
    }   
}

###
###
opfn <- paste(outdir, "Figure_S6_1_DEG_enrich.xlsx", sep="")
write.xlsx(plotDF, opfn)


###
### forest plots

fn <- paste(outdir, "Figure_S6_1_DEG_enrich.xlsx", sep="")
plotDF <- read.xlsx(fn)


col1 <- c("Bcell"="#4daf4a", "Monocyte"="#984ea3", "NKcell"="#aa4b56", "Tcell"="#ffaa00")
col2 <- c("LPS"="#fb9a99", "LPS+DEX"="#e31a1c", "PHA"="#a6cee3", "PHA+DEX"="#1f78b4")


plotDF <- plotDF%>%
    mutate(condition2=gsub("-", "+", conditions), MCls=gsub("_.*", "", conditions),
           contrast=gsub(".*_", "", conditions),
           log_odds=log(odds), log_lower=log(CI_lower), log_upper=log(CI_upper))

p <- ggplot(plotDF, aes(x=log_odds, y=condition2))+
   geom_errorbarh(aes(xmax=log_upper, xmin=log_lower, colour=MCls),
       size=0.5, height=0.2)+ 
   geom_point(aes(colour=MCls), shape=19, size=1.5)+
   scale_colour_manual(values=col1)+
   geom_vline(aes(xintercept=0), size=0.25, linetype="dashed")+ 
   xlab("log odds ratio")+xlim(-3, 1)+
   ###scale_y_discrete(labels=ylab2)+ 
   theme_bw()+
   theme(##plot.title=element_text(hjust=0.5, size=14),
         axis.title.y=element_blank(),
         axis.title.x=element_text(size=10),
         axis.text.x=element_text(size=10),
         axis.text.y=element_text(size=9),
         ## legend.position="none")
         legend.title=element_blank(),
         legend.text=element_text(size=9),
         legend.key.size=unit(0.4, "cm"))
         ## legend.position="none")
 
figfn <- paste(outdir, "Figure_S6_1_DEG_enrich.forest.png", sep="")
ggsave(figfn, p, device="png", width=550, height=480, units="px", dpi=120)    





###############################
### summary peak annotation 
##################################

fn <- paste(outdir, "TableS5_3_asthma-risk-genes_peak.txt.gz", sep="")
res <- fread(fn, header=T, sep="\t")

res <- res%>%mutate(comb=paste(MCls, contrast, sep="_"),
                    is_sig=ifelse(abs(estimate)>0.5&p.adjusted<0.1, 1, 0))
###
mat <- res%>%
    pivot_wider(id_cols=Gene, names_from=comb, values_from=is_sig, values_fill=0)%>%
    as.data.frame()

###
### 
fn <- paste(outdir, "TableS6_3_asthma-risk-genes_peak.txt.gz", sep="")
res <- fread(fn, header=T, sep="\t")

res <- res%>%mutate(comb=paste(MCls, contrast, sep="_"),
                    is_sig=ifelse(abs(estimate)>0.5&p.adjusted<0.1, 1, 0))



###
mat2 <- res%>%
    pivot_wider(id_cols=Gene, names_from=comb, values_from=is_sig, values_fill=0)%>%
    as.data.frame()

identical(colnames(mat), colnames(mat2))

mat_comb <- rbind(mat, mat2)

comb <- names(mat_comb)[-1]
comb <- comb[!grepl("^DC", comb)]
###
plotDF <- map_dfr(comb, function(ii){
   ##
   mat0 <- mat_comb[, c("Gene", ii)]
   names(mat0)[2] <- "is_sig"
   gene2 <- mat0%>%filter(is_sig>0)%>%pull(Gene)%>%unique()
   ## 
   plotDF0 <- data.frame(ny=length(gene2), conditions=ii)
   plotDF0
})


plotDF <- plotDF%>%filter(ny>0)%>%
    mutate(MCls=gsub("_.*", "", conditions), contrast=gsub(".*_", "", conditions),
           condition2=gsub("-", "+", conditions))


##
p0 <- ggplot(plotDF, aes(x=condition2, y=ny))+
    geom_bar(stat="identity", fill="#f768a1")+
    ## scale_fill_manual(values=col2,
    ##                   labels=c("Bcell"="B cell", "Monocyte"="Monocyte",
    ##                             "NKcell"="NK cell", "Tcell"="T cell"))+
    coord_flip()+
    ylab("#Genes in DARs")+    
    theme_bw()+
    theme(legend.title=element_blank(),
          axis.title.x=element_text(size=10),
          axis.title.y=element_blank(),
          axis.text=element_text(size=8))

figfn <- paste(outdir, "Figure_S6_3_DARs.bar.png", sep="")
ggsave(figfn, p0, width=380, height=450, units="px", dpi=120)



#######################
### plot 2 in DARs  ###
#######################


###
### aloft

fn <- paste(outdir, "TableS5_3_asthma-risk-genes_peak.txt.gz", sep="")
res <- fread(fn, header=T, sep="\t")

res <- res%>%mutate(comb=paste(MCls, contrast, sep="_"),
                    is_sig=ifelse(abs(estimate)>0.5&p.adjusted<0.1, 1, 0))

### GTEx
fn <- paste(outdir, "TableS6_3_asthma-risk-genes_peak.txt.gz", sep="")
res2 <- fread(fn, header=T, sep="\t")

res2 <- res2%>%mutate(comb=paste(MCls, contrast, sep="_"),
                    is_sig=ifelse(abs(estimate)>0.5&p.adjusted<0.1, 1, 0))

identical(colnames(res), colnames(res2))


###
### combine

res <- rbind(res, res2)


fn <- "/nfs/rprdata/julong/sc-atac/analyses.2021-02-05/2_Differential/1.3_DiffPeak.outs/3.0_DESeq_indi.results.rds"
resDP <- read_rds(fn)%>%
    mutate(comb=paste(MCls, contrast, sep="_"), is_sig=ifelse(abs(estimate)>0.5&p.adjusted<0.1, 1, 0))%>%
    as.data.frame()

comb2 <- sort(unique(res$comb))
comb2 <- comb2[!grepl("^DC", comb2)]

plotDF <- map_dfr(comb2, function(ii){
   ###
   gene_DAR <- res%>%filter(comb==ii, is_sig==1)%>%pull(peak)%>%unique() 
   DAR <- resDP%>%filter(comb==ii, is_sig==1)%>%pull(gene)%>%unique()
   df0 <- data.frame(conditions=ii, ngene_DAR=length(gene_DAR), nDAR=length(DAR))
   df0
})



plotDF2 <- plotDF%>%filter(ngene_DAR>0)%>%
    mutate(prop=(ngene_DAR/nDAR)*100, condition2=gsub("-", "+", conditions))


##
p0 <- ggplot(plotDF2, aes(x=condition2, y=prop))+
    geom_bar(stat="identity", fill="#f768a1")+
    ## scale_fill_manual(values=col2,
    ##                   labels=c("Bcell"="B cell", "Monocyte"="Monocyte",
    ##                             "NKcell"="NK cell", "Tcell"="T cell"))+
    coord_flip()+
    ylab("#Proportion of DARs (%)")+    
    theme_bw()+
    theme(legend.title=element_blank(),
          axis.title.x=element_text(size=10),
          axis.title.y=element_blank(),
          axis.text=element_text(size=8))

figfn <- paste(outdir, "Figure_S6_3_1_prop_DARs.bar.png", sep="")
ggsave(figfn, p0, width=380, height=450, units="px", dpi=120)


#############################
#### enrichment analysis  ###
#############################


fn <- paste(outdir, "TableS5_3_asthma-risk-genes_peak.txt.gz", sep="")
res <- fread(fn, header=T, sep="\t")

res <- res%>%mutate(comb=paste(MCls, contrast, sep="_"),
                    is_sig=ifelse(abs(estimate)>0.5&p.adjusted<0.1, 1, 0))

### GTEx
fn <- paste(outdir, "TableS6_3_asthma-risk-genes_peak.txt.gz", sep="")
res2 <- fread(fn, header=T, sep="\t")

res2 <- res2%>%mutate(comb=paste(MCls, contrast, sep="_"),
                    is_sig=ifelse(abs(estimate)>0.5&p.adjusted<0.1, 1, 0))

identical(colnames(res), colnames(res2))


###
### combine

res <- rbind(res, res2)
peakSel <- unique(res$peak)


fn <- "/nfs/rprdata/julong/sc-atac/analyses.2021-02-05/2_Differential/1.3_DiffPeak.outs/3.0_DESeq_indi.results.rds"
resDP <- read_rds(fn)%>%
    mutate(conditions=paste(MCls, contrast, sep="_"),
           is_sig=ifelse(abs(estimate)>0.5&p.adjusted<0.1, "is_DP", "not_DP"))%>%
    as.data.frame()

comb2 <- sort(unique(resDP$conditions))
comb2 <- comb2[!grepl("^DC", comb2)]

 
plotDF <- NULL

for( ii in comb2){
   ###
 
   tmp <- resDP%>%filter(conditions==ii)%>%
       mutate(is_twas_peak=ifelse(gene%in%peakSel, "is_asthma_peak", "not_asthma_peak"))

   DAR <- tmp%>%filter(is_sig=="is_DP")%>%pull(gene)%>%unique() 
   dat2 <- table(tmp$is_sig, tmp$is_twas_peak)
  
   ##
   if ( dat2[1,1]>0){ 
       fisher <- fisher.test(dat2) ##, alternative="greater") 
       fisher2 <- fisher.test(dat2, alternative="greater")
       
       df0 <- data.frame(conditions=ii, odds=fisher$estimate,
                      CI_lower=fisher$conf.int[1], CI_upper=fisher$conf.int[2],
                      pval=fisher2$p.value, olap=dat2[1,1], ndar=length(DAR)) 
       ##
       plotDF <- rbind(plotDF, df0)
    }   
}

###
###
opfn <- paste(outdir, "Figure_S6_3_DAR_enrich.xlsx", sep="")
write.xlsx(plotDF, opfn)    

###################
## Forest plots ###
####################

fn <- paste(outdir, "Figure_S6_3_DAR_enrich.xlsx", sep="")
plotDF <- read.xlsx(fn)


col1 <- c("Bcell"="#4daf4a", "Monocyte"="#984ea3", "NKcell"="#aa4b56", "Tcell"="#ffaa00")
col2 <- c("LPS"="#fb9a99", "LPS+DEX"="#e31a1c", "PHA"="#a6cee3", "PHA+DEX"="#1f78b4")


plotDF <- plotDF%>%
    mutate(condition2=gsub("-", "+", conditions), MCls=gsub("_.*", "", conditions),
           contrast=gsub(".*_", "", conditions),
           log_odds=log(odds), log_lower=log(CI_lower), log_upper=log(CI_upper))

range(c(plotDF$log_lower, plotDF$log_upper)) 

p <- ggplot(plotDF, aes(x=log_odds, y=condition2))+
   geom_errorbarh(aes(xmax=log_upper, xmin=log_lower, colour=MCls),
       size=0.5, height=0.2)+ 
   geom_point(aes(colour=MCls), shape=19, size=1.5)+
   scale_colour_manual(values=col1)+
   geom_vline(aes(xintercept=0), size=0.25, linetype="dashed")+ 
   xlab("log odds ratio")+xlim(-3, 3.5)+
   ###scale_y_discrete(labels=ylab2)+ 
   theme_bw()+
   theme(##plot.title=element_text(hjust=0.5, size=14),
         axis.title.y=element_blank(),
         axis.title.x=element_text(size=10),
         axis.text.x=element_text(size=10),
         axis.text.y=element_text(size=9),
         ## legend.position="none")
         legend.title=element_blank(),
         legend.text=element_text(size=9),
         legend.key.size=unit(0.4, "cm"))
         ## legend.position="none")
 
figfn <- paste(outdir, "Figure_S6_3_DAR_enrich.forest.png", sep="")
ggsave(figfn, p, device="png", width=550, height=480, units="px", dpi=120)    










#######################################
### summary motif annotation
#######################################

### aloft
fn <- paste(outdir, "TableS5_4_asthma-risk-genes_motif.txt.gz", sep="")
res <- fread(fn, header=T, sep="\t")

### gtex
fn <- paste(outdir, "TableS6_4_asthma-risk-genes_motif.txt.gz", sep="")
res2 <- fread(fn, header=T, sep="\t")

identical(colnames(res), colnames(res2))

## combine 
res <- rbind(res, res2)


###
### response motifs
contrast <- c("LPS", "LPS-DEX", "PHA", "PHA-DEX")
MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
combs <- paste(rep(MCls, each=4), rep(contrast, times=4), sep="_")
##
plotDF <- map_dfr(combs, function(ii){
   ##
   prefix <- "/nfs/rprdata/julong/sc-atac/analyses.2021-02-05/3_motif/2_motif.activities.outs/MotifList/" 
   fn <- paste(prefix, "1_comb_", ii, "_motif.txt", sep="")
   motif <- unique(read.table(fn)$V1)
   ### 
   res0 <- res%>%filter(motif_ID%in%motif)
   df0 <- data.frame(conditions=ii, ngene=length(unique(res0$Gene)), nmotif=length(unique(res0$motif_ID)))
   df0 
})    

plotDF <- plotDF%>%
    mutate(condition2=gsub("-", "+", conditions),
                          condition2=fct_reorder(condition2, ngene)) 

##
p0 <- ggplot(plotDF, aes(x=condition2, y=ngene))+
    geom_bar(stat="identity", fill="#f768a1")+
    ## scale_fill_manual(values=col2,
    ##                   labels=c("Bcell"="B cell", "Monocyte"="Monocyte",
    ##                             "NKcell"="NK cell", "Tcell"="T cell"))+
    coord_flip()+
    ylab("#Genes in Response motifs")+    
    theme_bw()+
    theme(legend.title=element_blank(),
          axis.title.x=element_text(size=10),
          axis.title.y=element_blank(),
          axis.text.x=element_text(size=8))

figfn <- paste(outdir, "Figure_S6_4_motifs.bar.png", sep="")
ggsave(figfn, p0, width=380, height=450, units="px", dpi=120)


#######################
### motif bar plots ###
#######################


### aloft
fn <- paste(outdir, "TableS5_4_asthma-risk-genes_motif.txt.gz", sep="")
res <- fread(fn, header=T, sep="\t")

### gtex
fn <- paste(outdir, "TableS6_4_asthma-risk-genes_motif.txt.gz", sep="")
res2 <- fread(fn, header=T, sep="\t")

identical(colnames(res), colnames(res2))

### combine aloft and gtex results
res <- rbind(res, res2)


### data for plot
summ <- res%>%group_by(Gene, motif_ID)%>%slice_min(pval, n=1)%>%ungroup()

plotDF <- summ%>%group_by(motif_name)%>%summarize(ngene=length(unique(Gene)), .groups="drop")%>%
    as.data.frame()

plotDF <- plotDF%>%arrange(desc(ngene))
plotDF2 <- plotDF%>%filter(ngene>5)
plotDF2$ny <- as.numeric(1:nrow(plotDF2))
plotDF2 <- plotDF2%>%mutate(motif_name2=fct_reorder(motif_name, ny, .desc=T))
 
 

motifSel <- plotDF2$motif_name


p0 <- ggplot(plotDF2, aes(x=motif_name2, y=ngene))+
    geom_bar(stat="identity", fill="#f768a1")+
    ## scale_fill_manual(values=col2,
    ##                   labels=c("Bcell"="B cell", "Monocyte"="Monocyte",
    ##                             "NKcell"="NK cell", "Tcell"="T cell"))+
    coord_flip()+
    ylab("#Genes in Response motifs")+    
    theme_bw()+
    theme(legend.title=element_blank(),
          axis.title.x=element_text(size=10),
          axis.title.y=element_blank(),
          axis.text.x=element_text(size=8),
          axis.text.y=element_text(size=8))

figfn <- paste(outdir, "Figure_S6_5_1_motifs.bar.png", sep="")
ggsave(figfn, p0, width=420, height=550, units="px", dpi=120)




################
### motif data 
################

dat <- res[,-c(1:3)]

dat2 <- res2[,-c(1:3)]

dat <- rbind(dat, dat2)

dat <- dat%>%distinct(conditions, motif_ID, .keep_all=T)
dat <- dat%>%group_by(conditions, motif_name)%>%slice_min(order_by=pval, n=1)%>%ungroup()


mat <- dat%>%
    pivot_wider(id_cols=motif_name, names_from=conditions, values_from=beta, values_fill=0)%>%
    column_to_rownames(var="motif_name")%>%as.matrix()



### significance
### response motifs


fn <- "/nfs/rprdata/julong/sc-atac/analyses.2021-02-05/3_motif/2_motif.activities.outs/3_motif.diff.results.rds"
resDiff <- read_rds(fn)%>%mutate(conditions=paste(MCls, contrast, sep="_"))
resDiff <- resDiff%>%dplyr::rename(motif_ID=gene, motif_name=motif)

contrast <- c("LPS", "LPS-DEX", "PHA", "PHA-DEX")
MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
combs <- paste(rep(MCls, each=4), rep(contrast, times=4), sep="_")
##

df_resp <- map_dfr(combs, function(ii){    
   ##
   prefix <- "/nfs/rprdata/julong/sc-atac/analyses.2021-02-05/3_motif/2_motif.activities.outs/MotifList/" 
   fn <- paste(prefix, "1_comb_", ii, "_motif.txt", sep="")
   motif <- unique(read.table(fn)$V1)
   motif2 <- resDiff%>%filter(conditions==ii, motif_ID%in%motif)%>%pull(motif_name)%>%unique() 
   df0 <- data.frame(conditions=ii, is_sig=1, motif_name=motif2)
   df0
})    

mat_sig <- df_resp%>%
    pivot_wider(id_cols=motif_name, names_from=conditions, values_from=is_sig, values_fill=0)%>%
    column_to_rownames(var="motif_name")%>%as.matrix()


mat <- mat[motifSel,]
mat_sig <- mat_sig[motifSel,]

identical(rownames(mat), rownames(mat_sig))
identical(colnames(mat), colnames(mat_sig))


mat2 <- mat*mat_sig



###
### Heatmap


###
### get colnames and re-order by treats
rn <- gsub("-", "+", colnames(mat2))
colnames(mat2) <- rn

x <- str_split(rn, "_", simplify=T)
cvt <- data.frame(comb=rn, MCls=x[,1], contrast=x[,2])
cvt <- cvt%>%arrange(contrast)

mat2 <- mat2[,cvt$comb]


###
### color for heatmap value
y <- as.vector(mat2)
## y0 <- y[abs(y)<2]
## quantile(abs(y), probs=0.99)

mybreak <- c(min(y,na.rm=T), seq(-1, 1, length.out=98), max(y,na.rm=T))

## quantile(abs(y), probs=c(0.9,0.95,0.99))
## range(y)

mycol <- colorRamp2(mybreak, colorRampPalette(rev(brewer.pal(n=7,name="RdBu")))(100))

 
###
### annotation columns
condition3 <- colnames(mat2)
col_ha <- HeatmapAnnotation(
   celltype=gsub("_.*", "", condition3),
   contrast=gsub(".*_", "", condition3),
   col=list(
      celltype=c("Bcell"="#4daf4a", "Monocyte"="#984ea3",
                 "NKcell"="#aa4b56", "Tcell"="#ffaa00"),
      contrast=c("LPS"="#fb9a99", "LPS+DEX"="#e31a1c",
                  "PHA"="#a6cee3", "PHA+DEX"="#1f78b4")),
   annotation_legend_param=list(celltype=list(labels_gp=gpar(fontsize=9), title_gp=gpar(fontsize=9),
      grid_width=grid::unit(0.38, "cm"), grid_height=grid::unit(0.45, "cm") ),
                                contrast=list(labels_gp=gpar(fontsize=9), title_gp=gpar(fontsize=9),
      grid_width=grid::unit(0.38, "cm"), grid_height=grid::unit(0.45, "cm"))),
   annotation_name_gp=gpar(fontsize=9),
   show_legend=c(T, T))
 
## x <- str_split(colnames(mat2), "_", simplify=T)
## df_col <- data.frame(celltype=x[,1], contrast=x[,2])
## col_ha <- HeatmapAnnotation(df=df_col, col=list(celltype=col1, contrast=col2),
##     annotation_name_gp=gpar(fontsize=10),
##     annotation_legend_param=list(
##   celltype=list(title_gp=gpar(fontsize=9), labels_gp=gpar(fontsize=9), title="celltype",
##                 grid_width=unit(0.45, "cm"), grid_height=unit(0.5, "cm")),
##   contrast=list(title_gp=gpar(fontsize=9), labels_gp=gpar(fontsize=9), title="contrast",
##                 grid_width=unit(0.45, "cm"), grid_height=unit(0.5, "cm"))))

###
### 111 motifs
mat2[mat2==0] <- NA
p0 <- Heatmap(mat2, col=mycol, na_col="white",
   ###rect_gp=gpar(col="black", lwd=1),           
   cluster_rows=F,  cluster_columns=F,
   show_row_dend=F, show_column_dend=F, 
   show_row_names=T, row_names_gp=gpar(fontsize=7), row_names_side="left",
   show_column_names=T, column_names_gp=gpar(fontsize=7),
   column_names_rot=-45,   
   top_annotation=col_ha,
   ###show_heatmap_legend=F,
   heatmap_legend_param=list(title="Diff motif",
      title_gp=gpar(fontsize=9),
      at=seq(-1, 1, by=0.5), 
      labels_gp=gpar(fontsize=9),
      grid_width=grid::unit(0.38, "cm"),
      legend_height=grid::unit(5, "cm")))
 
 
###
figfn <- paste(outdir, "Figure_S6_5_2_condition_motif.heatmap.png", sep="")
png(figfn, width=520, height=600,res=120)
set.seed(0)
p0 <- draw(p0)
dev.off()
 




###
### alter

fn <- "./5_pub.outs/2_supp_tables/TableS5_2_asthma-risk-genes_in_response_ALOFT.txt.gz"
res <- fread(fn, header=T, sep="\t", data.table=F)

###
fn2 <- "./5_pub.outs/2_supp_tables/TableS6_2_asthma-risk-genes_in_response.txt.gz"
res2 <- fread(fn2, header=T, sep="\t", data.table=F)

## 
res_comb <- rbind(res, res2)

### 
score_diff <- sapply(1:nrow(res_comb), function(i){    
   ###
   motif <- unlist(str_split(res_comb$motifAnnot[i], ";"))
   motif <- str_split(motif, "_", simplify=T)
   motif2 <- data.frame(motif_ID=motif[,1], score_ref=as.numeric(motif[,2]), score_alt=as.numeric(motif[,3]))%>%
       mutate(diff=abs(score_ref-score_alt))
   s0 <- max(motif2$diff) 
   s0
})    

res_comb <- cbind(res_comb, score_diff)

th0 <- 3*log2(exp(1))    
res_sig <- res_comb%>%filter(score_diff>th0)




















###
###
fn <- paste(outdir, "TableS5_2_asthma-risk-genes_in_response.txt.gz", sep="")
res <- fread(fn, header=T, sep="\t", data.table=F)

## aloft
fn <- "/nfs/rprdata/julong/sc-atac/eqtl_analysis_plots/2_tables/TableS4_3_eGene_ALOFT_dap.txt.gz"
dap_aloft <- read.table(fn, header=T)
egene_aloft <- unique(dap_aloft$gene)

res <- res%>%
    mutate(is_egene_aloft=ifelse(Gene%in%dap_aloft$gene, 1, 0),
           is_egene_sc=ifelse(Gene%in%dap_sc$gene, 1, 0))

## ### scaip
## fn <- "/nfs/rprdata/julong/sc-atac/eqtl_analysis_plots/2_tables/TableS4_2_eGene_scaip_dap.txt.gz"
## dap_sc <- fread(fn, header=T)%>%dplyr::select(gene, lfdr, conditions)
## egene_sc <- unique(dap_sc$gene)

## dap_sc$is_sig <- 1
## mat_sc <- dap_sc%>%
##     pivot_wider(id_cols=gene, names_from=conditions, values_from=is_sig, values_fill=0)

## mat_sc$counts <- rowSums(mat_sc[,2:21])

## res <- res

fn <- "/nfs/rprdata/julong/SCAIP/analyses/SCAIP-B1-6_2020.03.23/6_DEG.CelltypeNew_output/Filter2/2_meta.rds"
res_diff <- read_rds(fn)%>%
    mutate(zscore_diff=beta/stderr, conditions=paste(MCls, contrast, sep="_"),
           is_sig=ifelse(qval<0.01, 1, 0),
           zscore_diff2=zscore_diff*is_sig)
 
sigs <- res_diff%>%filter(is_sig==1)%>%pull(gene)%>%unique()
  
mat_diff <- res_diff%>%
    pivot_wider(id_cols=gene, names_from=conditions, values_from=zscore_diff2, values_fill=0) 


###
###
 
res2 <- res%>%dplyr::select(Gene, symbol, zscore_twas, PIP)%>%
    left_join(mat_diff, by=c("Gene"="gene"))
## res2 <- res2%>%filter(Gene%in%sig)

opfn <- paste(outdir, "TableS5_6_asthma-risk-genes_DEG.xlsx", sep="")
write.xlsx(res2, file=opfn)



#################
#### plots 
#################


outdir2 <- "./5_pub.outs/3_example_plots/tmp_plots/"
if (!file.exists(outdir2)) dir.create(outdir2, showWarnings=F, recursive=T)


fn <- "./5_pub.outs/3_example_plots/TableS5_6_asthma-risk-genes_DEG.xlsx"
res <- read.xlsx(fn)

plotDF <- res%>%
    pivot_longer(cols=!(Gene|symbol|zscore_twas|PIP), names_to="conditions", values_to="zscore_diff")

plotDF <- plotDF%>%
    mutate(MCls=gsub("_.*", "", conditions), contrast=gsub(".*_", "", conditions))%>%
    drop_na(zscore_diff)

sigs <- plotDF%>%filter(zscore_diff!=0, PIP>0.1)%>%pull(Gene)%>%unique()

####
###
plotDF2 <- plotDF%>%dplyr::filter(grepl("DEX", contrast), zscore_diff!=0, PIP>0.1)
p2 <- ggplot(plotDF2, aes(x=zscore_diff, y=zscore_twas, color=factor(MCls)))+
   geom_point(shape=24)+
   geom_text_repel(data=plotDF2,
      mapping=aes(x=zscore_diff, y=zscore_twas, label=symbol, color=factor(MCls)), 
      fontface="italic", size=2.5, max.overlaps=Inf,
      min.segment.length=0, segment.curvature=-0.1, segment.angle=20)+  
      scale_color_manual(values=c("Bcell"="#4daf4a", "Monocyte"="#984ea3",
         "NKcell"="#aa4b56", "Tcell"="#ffaa00"), guide="none")+
      xlab(bquote(~italic(Z)~"-score of DE"))+
      ylab(bquote(~italic(Z)~"-score of asthma TWAS"))+
      xlim(-15, 36)+ylim(-11, 11)+
      geom_hline(yintercept=0, linetype="dashed", color="grey30")+
      geom_vline(xintercept=0, linetype="dashed", color="grey30")+
      facet_grid(MCls~contrast, 
                 labeller=labeller(MCls=c("Bcell"="B cell", "Monocyte"="Monocyte",
                                          "NKcell"="NK cell", "Tcell"="T cell"),
                                  contrast=c("PHA"="LPS", "PHA"="PHA",
                                             "LPS-DEX"="LPS+DEX", "PHA-DEX"="PHA+DEX")))+
      theme_bw()+
      theme(## legend.position=c(0.25, 0.85),
            ## legend.title=element_blank(),
            ## legend.background=element_blank(),
            ## legend.box.background=element_blank(),
            ## legend.key=element_blank(),
            ## plot.title=element_text(hjust=0.5),
            strip.text=element_text(size=12),
            axis.title=element_text(size=10))
###
figfn <- paste(outdir2, "Figure1.2_twas_DEG_DEX.png", sep="")
ggsave(figfn, p2, width=520, height=720, units="px", dpi=120)



###
### Not DEX 
plotDF2 <- plotDF%>%dplyr::filter(!grepl("DEX", contrast), zscore_diff!=0, PIP>0.1)
p0 <- ggplot(plotDF2, aes(x=zscore_diff, y=zscore_twas, color=factor(MCls)))+
   geom_point(shape=24)+
   geom_text_repel(data=plotDF2,
      mapping=aes(x=zscore_diff, y=zscore_twas, label=symbol, color=factor(MCls)), 
      fontface="italic", size=2.5, max.overlaps=Inf,
      min.segment.length=0, segment.curvature=-0.1, segment.angle=20)+  
      scale_color_manual(values=c("Bcell"="#4daf4a", "Monocyte"="#984ea3",
         "NKcell"="#aa4b56", "Tcell"="#ffaa00"), guide="none")+
      xlab(bquote(~italic(Z)~"-score of DE"))+
      ylab(bquote(~italic(Z)~"-score of asthma TWAS"))+
      xlim(-15, 36)+ylim(-11, 11)+
      geom_hline(yintercept=0, linetype="dashed", color="grey30")+
      geom_vline(xintercept=0, linetype="dashed", color="grey30")+
      facet_grid(MCls~contrast, 
                 labeller=labeller(MCls=c("Bcell"="B cell", "Monocyte"="Monocyte",
                                          "NKcell"="NK cell", "Tcell"="T cell"),
                                  contrast=c("LPS"="LPS", "PHA"="PHA",
                                             "LPS-DEX"="LPS+DEX", "PHA-DEX"="PHA+DEX")))+
      theme_bw()+
      theme(## legend.position=c(0.25, 0.85),
            ## legend.title=element_blank(),
            ## legend.background=element_blank(),
            ## legend.box.background=element_blank(),
            ## legend.key=element_blank(),
            ## plot.title=element_text(hjust=0.5),
            strip.text=element_text(size=12),
            axis.title=element_text(size=10))
###
figfn <- paste(outdir2, "Figure1.1_twas_DEG.png", sep="")
ggsave(figfn, p0, width=520, height=720, units="px", dpi=120)

sigs <- plotDF%>%filter(zscore_diff!=0, PIP>0.1)%>%pull(Gene)%>%unique()



###
###





fn <- "./4_INTACT.outs/Asthma_Bothsex_inv_var_meta_GBMI_052021_nbbkgt1_gtex_topPIP_union_combinfor.txt"
res <- fread(fn, header=T, data.table=F, sep="\t")
x <- str_split(res$genetic_variant, ":", simplify=T)
res$pos <- as.integer(x[,2])


##
fn <- "/nfs/rprdata/julong/sc-atac/twas_analysis_2022-10-16/analysis2/2_gwas_prepare/SCAIP_final_bed.gz"
bed <- fread(fn, header=T, sep="\t")

query_snp <- c("rs7216389", "rs8076131", "rs2305480")

bed2 <- bed%>%filter(rs%in%query_snp)


ss <- 100000
for ( i in 1:nrow(bed2)){
    ##
    chr_i <- as.integer(gsub("_.*", "", bed2$chr_pos_grch38[i]))
    pos_i <- as.integer(gsub(".*_", "", bed2$chr_pos_grch38[i]))
 
    res2 <- res%>%filter(chromosome==chr_i, pos>(pos_i-ss), pos<(pos_i+ss), FDR_intact<0.1)%>%
        mutate(diff=abs(pos-pos_i))
    
    k <- 1
    ens0 <- res2$Gene[k]
    snp_id <- bed2$chr_pos_grch38
    
    fn <- "/nfs/rprdata/julong/sc-atac/twas_analysis_2022-10-16/Whole_Blood_GTEx_v8_results/Whole_Blood_allSNPs_union.txt.gz"
    dap <- fread(fn, header=F, data.table=F)
    names(dap) <- c("gene", "genetic_variant", "PIP")
    dap <- dap%>%mutate(gene=gsub("\\..*", "", gene))
    dap2 <- dap%>%filter(gene==ens0)
    x <- str_split(dap2$genetic_variant, "_", simplify=T)
    dap2$chr_pos_grch38 <- paste(gsub("chr", "", x[,1]), x[,2], sep="_") 
    
    dap2 <-dap2%>%filter(chr_pos_grch38%in%snp_id)
 







    
###################
### cluster PIP ###
###################
 
## sub_i <- c(7, 11, 12, 14, 20, 21, 25)

## res3 <- res2%>%dplyr::filter(symbol%in%c("GSDMB", "ORMDL3", "IKZF3", "IL4", "FADS2"))

i <- 5   
ens <- res3$gene[i]
symbol <- res3$symbol[i]
symbol2 <- symbol    
chr_pos <- res3$chr_pos_grch37[i]
chr_i <- as.numeric(gsub("_.*", "", chr_pos))
pos_i <- as.numeric(gsub(".*_", "",  chr_pos))
snp_i <- res3$genetic_variant[i]

cat(ens, symbol2, "\n")




prefix <- "/nfs/rprdata/julong/sc-atac/genetic_analysis_ALOFT/DAP-G/dap-g_outs/dap-g_combineNew/Union/"
fn <- paste(prefix, ens, ".SNP.out", sep="")
tmp <- read.table(fn)
tmp2 <- tmp%>%filter(V5==1)
###
sum(tmp2$V3)
nrow(tmp2)


###
###
df <- res3[sub_i, c("gene", "symbol")]
fn <- "./enloc_analysis/ALOFT_intact.txt"
x <- read.table(fn, header=T)
x2 <- x%>%filter(Gene%in%df$gene)%>%left_join(df, by=c("Gene"="gene"))


#####################
### response eQTL ###
#####################

load("/wsu/home/groups/piquelab/SCAIP/SCAIP-genetic/mashr_eQTL/mashr-reQTLs_union-unshared-magnitude2-mlfsr0.1.Rd")
names(uum) <- gsub("-EtOH", "", gsub(":.*", "", names(uum)))

reqtl <- map_dfr(names(uum), function(ii){
   ###
   gene_SNP <- uum[[ii]] 
   tmp <- data.frame(gene=gsub("\\..*", "", gsub("_.*", "", gene_SNP)),
                     genetic_variant=gsub(".*_", "", gene_SNP), conditions=ii)
   tmp
})

## MICB
i <- 2
chr_pos <- res3$chr_pos_grch37[i]
chri <- as.numeric(gsub("_.*", "", chr_pos))
posi <- as.integer(gsub(".*_", "", chr_pos))
symbol_i <- res3$symbol[i]
SNP_i <- res3$genetic_variant[i]    


######################
### colocalization ###
######################

## res3 <- res2%>%filter(peaking_d==3, FDR_intact<0.1)

i <- 1
ens <- res3$gene[i]
symbol <- res3$symbol[i]
symbol2 <- symbol    
chr_pos <- res3$chr_pos_grch37[i]
chr_i <- as.numeric(gsub("_.*", "", chr_pos))
pos_i <- as.numeric(gsub(".*_", "",  chr_pos))
snp_i <- res3$genetic_variant[i]


fn <- "./enloc_analysis/enloc.gene.out"
enloc_gene <- read.table(fn, header=T)%>%arrange(desc(GLCP))

gene2 <- x%>%filter(GLCP>0.5)%>%pull(Gene)
 

enloc_snp <- read.table("./enloc_analysis/enloc.snp.out", header=T)
x <- enloc_snp %>%filter(SNP==snp_i)

###
intact <- read.table("./enloc_analysis/ALOFT_intact.txt", header=T)
                      
## %>%
##     mutate(gene=gsub(":.*", "", Signal))
## x%>%filter(gene==ens)


gwas <- fread("./enloc_analysis/Asthma_gwas_grch37.pip.gz", header=F)
