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

options(scipen=16)

rm(list=ls())


outdir <- "./5_pub.outs/2_supp_tables/"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)




#######################
### supp table 
########################
trait <- sort(read.table("traits.txt")$V1)[3]

fn <- paste("./4_INTACT.outs/", trait, "_aloft_topPIP_union_combinfor.txt", sep="")
res <- read.table(fn, header=T, sep="\t")

res2 <- res%>%
    dplyr::select(Gene, symbol, genetic_variant, chromosome=chr, chr_pos_grch38, PIP, pval_eqtl,
    peaking_d, pval_twas=pval_gwas, zscore_twas=zscore_gwas, FDR_twas, GLCP, PCG, FDR_intact=FDR)

res2 <- res2%>%arrange(FDR_intact)

###
opfn <- gzfile(paste(outdir, "TableS5_1_asthma-risk-genes_ALOFT.txt.gz", sep=""))
write.table(res2, file=opfn, quote=F, row.names=F, sep="\t")


####
fn <- paste(outdir, "TableS5_1_asthma-risk-genes_ALOFT.txt.gz", sep="")
res <- fread(fn, header=T, data.table=F)

x <- str_split(res$genetic_variant, ":", simplify=T)

res2 <- data.frame(res[,1:5], chr_pos_grch37=paste(x[,1], x[,2], sep="_"), res[,6:14])
opfn <- gzfile(paste(outdir, "TableS5_1_asthma-risk-genes_ALOFT.txt.gz", sep=""))
write.table(res2, file=opfn, quote=F, row.names=F, sep="\t")




####################################
### annotation genetic variants 
#####################################


fn <- paste(outdir, "TableS5_1_asthma-risk-genes_ALOFT.txt.gz", sep="")
res <- fread(fn, header=T, data.table=F)

res2 <- res%>%filter(FDR_intact<0.1, peaking_d==3, PIP>0.1) 

###
bed <- str_split(res2$genetic_variant, ":", simplify=T)
res2 <- res2%>%mutate(pos=as.integer(bed[,2]))




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
   chri <- res2$chromosome[i]
   posi <- res2$pos[i] 
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


## x <- snpmotif2%>%filter(abs(score_alt-score_ref)>log2(20))
## snp_sig <- unique(x$variants)
## ##
## fn <- "/nfs/rprdata/julong/sc-atac/genetic_analysis_ALOFT/DAP-G/5_summary.outs/aloft_allSNPs_union.txt.gz"
## dap <- fread(fn, header=F)
## snp2 <- dap%>%filter(V3>0.9)%>%pull(V2)%>%unique()

## snp2_sig <- intersect(snp_sig, snp2)

###
annoMotif <- sapply(1:nrow(res2), function(i){    
  ##   
  chr_pos_i <- paste(res2$chromosome[i], res2$pos[i], sep="_") 
  ##  
  motifSel <- snpmotif2%>%filter(chr_pos==chr_pos_i)%>%
      mutate(motif_score=paste(motifs, round(score_ref, 3), round(score_alt, 3), sep="_"))%>%
      pull(motif_score)  
  ###   
  motif0 <- paste(motifSel, collapse=";")
  motif0  
})    
 
res2$motifAnnot <- as.character(annoMotif)

### txt.gz format
opfn <- gzfile(paste(outdir, "TableS5_2_asthma-risk-genes_in_response.txt.gz", sep=""))
write.table(res2, file=opfn, quote=F, row.names=F, sep="\t")

### save xlsx file 
opfn <- paste(outdir, "TableS5_2_asthma-risk-genes_in_response.xlsx", sep="")
write.xlsx(res2, file=opfn)


#################
### ase 
#################

fn <- paste(outdir, "TableS5_2_asthma-risk-genes_in_response.txt.gz", sep="")
res <- fread(fn, header=T, sep="\t")
res <- res%>%mutate(chr_pos=paste(chromosome, pos, sep="_"))

res <- res[,1:18]

###
fn2 <- "/nfs/rprdata/julong/sc-atac/demux.2021-01-23/asequant/asefinal.txt.gz"
ase <- fread(fn2, header=T, data.table=F)  
ase <- ase%>%mutate(chr_pos=paste(chr, pos, sep="_"))


###
###
ase_summ <- map_dfr(1:nrow(res), function(i){
   ##
   pos0 <- res$chr_pos[i]
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

res2 <- cbind(res, ase_summ)

### txt.gz format
opfn <- gzfile(paste(outdir, "TableS5_2_asthma-risk-genes_in_response.txt.gz", sep=""))
write.table(res2, file=opfn, quote=F, row.names=F, sep="\t")


### save xlsx file 
opfn <- paste(outdir, "TableS5_2_asthma-risk-genes_in_response.xlsx", sep="")
write.xlsx(res2, file=opfn)


###
### reorder column for final manuscript 
fn <- paste(outdir, "TableS5_2_asthma-risk-genes_in_response.txt.gz", sep="")
res <- read.table(fn, header=T, sep="\t")

x <- str_split(res$genetic_variant, ":", simplify=T)

res2 <- data.frame(res[,1:5], chr_pos_grch37=paste(x[,1], x[,2], sep="_"), res[, c(6:14, 16:17, 19:22)])
opfn <- gzfile(paste(outdir, "TableS5_2_asthma-risk-genes_in_response_ALOFT.txt.gz", sep=""))
write.table(res2, file=opfn, quote=F, row.names=F, sep="\t")


##########################
### differential
##########################

 
fn <- paste(outdir, "TableS5_2_asthma-risk-genes_in_response_ALOFT.txt.gz", sep="")
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
opfn <- gzfile(paste(outdir, "TableS5_3_asthma-risk-genes_peak.txt.gz", sep=""))
write.table(DF, file=opfn, quote=F, row.names=F, sep="\t")



####
### save xlsx fil
fn <- paste(outdir, "TableS5_3_asthma-risk-genes_peak.txt.gz", sep="")
res <- read.table(fn, header=T, sep="\t")

opfn <- paste(outdir, "TableS5_3_asthma-risk-genes_peak.xlsx", sep="")
write.xlsx(res, file=opfn)




    
###########################
### differential motif
###########################

fn <- paste(outdir, "TableS5_2_asthma-risk-genes_in_response.txt.gz", sep="")
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

opfn <- gzfile(paste(outdir, "TableS5_4_asthma-risk-genes_motif.txt.gz", sep=""))
write.table(DF, file=opfn, quote=F, row.names=F, sep="\t")


###
### save xlsx
fn <- paste(outdir, "TableS5_4_asthma-risk-genes_motif.txt.gz", sep="")
res <- read.table(fn, header=T, sep="\t")

opfn <- paste(outdir, "TableS5_4_asthma-risk-genes_motif.xlsx", sep="")
write.xlsx(res, file=opfn)






##################################
#### risk genes from GTEx
##################################












##############################
### summary ase results
###############################

## fn <- gzfile(paste(outdir, "TableS5_2_asthma-risk-genes_in_response.txt.gz", sep=""))
## res <- read.table(fn, header=T, sep="\t")
 
## fn <- "/nfs/rprdata/julong/sc-atac/demux.2021-01-23/asequant/asefinalsigAlt.txt.gz"
## ase <- fread(fn, header=T, data.table=F)  

## ase2 <- ase%>%filter(id%in%res$genetic_variant) 

## ase2 <- map_dfr(1:nrow(ase2), function(i){
##    ##
##    res0 <- res%>%dplyr::filter(genetic_variant==ase2$id[i])
##    tmp0 <- cbind(ase2[i,], res0[1,])
##    tmp0
## })

## opfn <- paste(outdir, "TableS5_5_asthma-risk-genes_ase.xlsx", sep="")
## write.xlsx(ase2, opfn)



########################
### extract results 
##########################

## outdir2 <- "./5_pub.outs/3_example_plots/Example2_table/"
## if (!file.exists(outdir2)) dir.create(outdir2, showWarnings=F, recursive=T)


## fn <- paste(outdir, "TableS5_2_asthma-risk-genes_in_response.txt.gz", sep="")
## res <- fread(fn, header=T, sep="\t")
## res2 <- res%>%filter(PIP>0.1, FDR_twas<0.1)%>%arrange(pval_twas)

## ### DP
## fn <- paste(outdir, "TableS5_3_asthma-risk-genes_peak.txt.gz", sep="")
## resDP <- fread(fn, header=T, sep="\t")


## ### motifs
## fn <- paste(outdir, "TableS5_4_asthma-risk-genes_motif.txt.gz", sep="")
## resMotif <- fread(fn, header=T, sep="\t")

## ###
## ###
## for (i in c(1:10, 50)){

## ###    
## symbol0 <- res2$symbol[i]
## ens <- res2$Gene[i]

## ### DARs
## resDP2 <- resDP%>%filter(symbol==symbol0)%>%
##     mutate(conditions=paste(MCls, contrast, sep="_"),
##            baseMean=round(baseMean, 3), estimate=round(estimate, 3),
##            p.value=round(p.value, 3), p.adjusted=round(p.adjusted, 3),
##            is_sig=ifelse(p.adjusted<0.1&abs(estimate)>0.5, 1, 0))%>%
##     dplyr::select(conditions, peak, estimate, pval=p.value, FDR=p.adjusted, is_sig)%>%
##     arrange(desc(abs(estimate)))

## opfn <- paste(outdir2, "Table_", i, "_", ens, "_", symbol0, "_DAR.xlsx", sep="")
## write.xlsx(resDP2, file=opfn)

    
## ### Motifs
## resMotif2 <- resMotif%>%filter(symbol==symbol0)%>%
##     mutate(conditions=paste(MCls, contrast, sep="_"),
##            beta=round(beta, 3), pval=round(pval, 3), qval=round(qval,3))%>%
##    dplyr::select(conditions, motif_ID, motif_name, beta, pval, qval)%>%
##    arrange(desc(abs(beta)))

## opfn2 <- paste(outdir2, "Table_", i, "_", ens, "_", symbol0, "_motif.xlsx", sep="")
## write.xlsx(resMotif2, file=opfn2)

## cat(symbol0, "\n")

## }    


#############################################
### summary the results of risk genes  
#############################################

  
fn <- paste(outdir, "TableS5_1_asthma-risk-genes_ALOFT.txt.gz", sep="")
res <- fread(fn, header=T, sep="\t")
res2 <- res%>%filter(FDR_intact<0.1)


## ### eGenes
## fn <- "/nfs/rprdata/julong/sc-atac/eqtl_analysis_plots/2_tables/TableS4_2_eGene_scaip_dap.txt.gz"
## dap <- fread(fn, header=T)
## dap <- dap%>%group_by(conditions, gene)%>%distinct(gene)

## plotDF <- map_dfr(sort(unique(dap$conditions)), function(ii){
##     ##
##     egene <- dap%>%filter(conditions==ii)%>%pull(gene)
##     olap <- intersect(egene, res2$Gene)
##     df0 <- data.frame(conditions=ii, olap=length(olap))
##     df0
## })    

## plotDF <- plotDF%>%filter(olap>0)%>%
##    mutate(condition2=gsub("-", "+", conditions),
##           condition2=fct_reorder(condition2, olap))

## ##
## p0 <- ggplot(plotDF, aes(x=condition2, y=olap))+
##     geom_bar(stat="identity", fill="#f768a1")+
##     coord_flip()+
##     ylab("#eGenes")+    
##     theme_bw()+
##     theme(legend.title=element_blank(),
##           axis.title.x=element_text(size=10),
##           axis.title.y=element_blank(),
##           axis.text=element_text(size=8))

## figfn <- paste(outdir, "Figure_S5_1_eGenes.bar.png", sep="")
## ggsave(figfn, p0, width=380, height=450, units="px", dpi=120)



### DEGs
fn <- "/nfs/rprdata/julong/SCAIP/analyses/SCAIP-B1-6_2020.03.23/6_DEG.CelltypeNew_output/Filter2/2_meta.rds"
res_diff <- read_rds(fn)%>%filter(abs(beta)>0.5&qval<0.1)%>%
    mutate(conditions=paste(MCls, contrast, sep="_"))%>%
    dplyr::select(conditions, gene)


plotDF <- map_dfr(sort(unique(res_diff$conditions)), function(ii){
    ##
    DEG <- res_diff%>%filter(conditions==ii)%>%pull(gene)
    olap <- intersect(DEG, res2$Gene)
    df0 <- data.frame(conditions=ii, olap=length(olap))
    df0
})    

plotDF <- plotDF%>%filter(olap>0)%>%
   mutate(condition2=gsub("-", "+", conditions),
          condition2=fct_reorder(condition2, olap))
 
##
p0 <- ggplot(plotDF, aes(x=condition2, y=olap))+
    geom_bar(stat="identity", fill="#f768a1")+
    coord_flip()+
    ylab("#DEGs")+    
    theme_bw()+
    theme(legend.title=element_blank(),
          axis.title.x=element_text(size=10),
          axis.title.y=element_blank(),
          axis.text=element_text(size=8))

figfn <- paste(outdir, "Figure_S5_1_DEGs.bar.png", sep="")
ggsave(figfn, p0, width=380, height=450, units="px", dpi=120)





###
### summary peak annotation 
fn <- paste(outdir, "TableS5_3_asthma-risk-genes_peak.txt.gz", sep="")
res <- fread(fn, header=T, sep="\t")

res <- res%>%mutate(comb=paste(MCls, contrast, sep="_"),
                    is_sig=ifelse(abs(estimate)>0.5&p.adjusted<0.1, 1, 0))

###
mat <- res%>%pivot_wider(id_cols=Gene, names_from=comb, values_from=is_sig, values_fill=0)

ngenes <- colSums(mat[,-1], na.rm=T)

plotDF <- data.frame("ny"=ngenes, conditions=names(ngenes))%>%filter(ny>0)%>%
    mutate(MCls=gsub("_.*", "", conditions), contrast=gsub(".*_", "", conditions),
           condition2=gsub("-", "+", conditions),
           condition2=fct_reorder(condition2, ny))

### setting color
col1 <- c("LPS"="#fb9a99", "LPS+DEX"="#e31a1c",
   "PHA"="#a6cee3", "PHA+DEX"="#1f78b4")
col2 <- c("Bcell"="#4daf4a", "Monocyte"="#984ea3",
   "NKcell"="#aa4b56", "Tcell"="#ffaa00")


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

figfn <- paste(outdir, "Figure_S5_3.bar.png", sep="")
ggsave(figfn, p0, width=380, height=450, units="px", dpi=120)


###
### summary motif annotation
 
fn <- paste(outdir, "TableS5_4_asthma-risk-genes_motif.txt.gz", sep="")
res <- fread(fn, header=T, sep="\t")
 
summ <- res%>%group_by(Gene, motif_name)%>%slice_min(pval, n=1)%>%ungroup()
summ <- summ%>%group_by(Gene, motif_name)%>%distinct(motif_name, .keep_all=T)%>%ungroup()

summ$is_in <- 1

mat <- summ%>%pivot_wider(id_cols=Gene, names_from=motif_name, values_from=is_in, values_fill=0)

ngenes <- colSums(mat[,-1])

plotDF <- data.frame(ny=ngenes, motif_name=names(ngenes))%>%arrange(desc(ny))
plotDF <- plotDF%>%mutate(motif_name2=fct_reorder(motif_name, ny))

##
p0 <- ggplot(plotDF%>%filter(ny>3), aes(x=motif_name2, y=ny))+
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
          axis.text=element_text(size=8))

figfn <- paste(outdir, "Figure_S5_4.bar.png", sep="")
ggsave(figfn, p0, width=380, height=450, units="px", dpi=120)


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
