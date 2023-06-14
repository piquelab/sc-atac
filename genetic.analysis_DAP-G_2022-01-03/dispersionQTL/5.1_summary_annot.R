###
library(Matrix)
library(tidyverse)
library(data.table)
library(clusterProfiler)
library(org.Hs.eg.db)
library(annotables)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(openxlsx)

###
### summary and annotation fine-mapping vQTL 
### By JW, May-16-2023

##################################
### find the annotation of SNP ###
##################################

###
### dap-g results

conditions <- read.table("datasets.txt")$V1

condition2 <- conditions[grepl("DEX", conditions)]


vgene <- lapply(conditions, function(ii){    
###    
condition <- ii
fn <- paste("./5_summary.outs/dap-g_pct_0.02_union/summary/", condition, "_gene_FDR.txt", sep="")
dap <- read.table(fn, header=T)%>%filter(FDR<0.1)
dap$gene
})
vgene <- unique(do.call(c, vgene))



res <- lapply(conditions, function(ii){
###    
condition <- ii
fn <- paste("./5_summary.outs/dap-g_pct_0.02_union/summary/", condition, "_gene_FDR.txt", sep="")
dap <- read.table(fn, header=T)%>%filter(FDR<0.1)    

cat(ii, "\n")
    
##
res_dap <- NULL    
if ( nrow(dap)>0){
   ##
   for ( i in 1:nrow(dap)){
       
      ens <- dap$gene[i] 
      fn <- paste("./dap-g_outs/dap-g_pct_0.02_union/", condition, "/", ens, ".SNP.out", sep="")
      res <- read.table(fn, header=F)%>%dplyr::filter(V3>0.1)
      snp1 <- res$V2
      pip1 <- as.numeric(res$V3)
      nsnp <- nrow(res) 
      ###
      tmp <- data.frame(conditions=rep(condition, nsnp), gene=rep(ens, nsnp), SNPs=snp1, pip=pip1)
      res_dap <- rbind(res_dap, tmp)
   }
}
###
res_dap
})

res <- do.call(rbind, res)
rownames(res) <- NULL

res <- res%>%mutate(MCls=gsub("_.*", "", conditions), treats=gsub(".*_", "", conditions))


### add annotation category
MCls <- sort(unique(res$MCls))
res <- map_dfr(MCls, function(oneMCl){
   ### SNP annotation
   cat(oneMCl, "\n") 
   prefix <- "/nfs/rprdata/julong/sc-atac/genetic.analysis_torus_2021-10-14/SNPannotation/"
   fn <- paste(prefix, "4_SNPAnnot.outs/pct_0.02/3_", oneMCl, "_union_torus.annot.gz", sep="")
   anno <- fread(fn, header=T, data.table=F)

   ###
   res2 <- res%>%dplyr::filter(MCls==oneMCl)
   res2 <- res2%>%left_join(anno, by=c("SNPs"="SNP")) 
   res2
})


###
###
geneID <- bitr(res$gene, fromType="ENSEMBL", toType="SYMBOL", OrgDb=org.Hs.eg.db)

res2 <- res%>%left_join(geneID, by=c("gene"="ENSEMBL"))

res2 <- res2%>%mutate(snpinfor=gsub(";.*", "", SNPs),
                      chr=gsub(":.*", "", snpinfor),
                      pos=gsub(":.*", "",  gsub("^\\d{1,2}:", "", snpinfor)),
                      chr_pos=paste(chr, pos, sep="_"))

 
opfn <- "./5_summary.outs/dap-g_pct_0.02_union/vQTL_annoted_summary.xlsx"
write.xlsx(res2, file=opfn, overwrite=T)




###
### combine information-add vQTL information 

fn <- "./5_summary.outs/dap-g_pct_0.02_union/vQTL_annoted_summary.xlsx"
res <- read.xlsx(fn)
res <- res%>%mutate(gene_SNP=paste(gene, SNPs, sep="_"))
comb <- unique(res$conditions)

###
resNew <- map_dfr(comb, function(ii){    
   ##
   res2 <- res%>%filter(conditions==ii)   
   cat(ii, "\n") 
    
   fn <- paste("./eQTL_results/", ii, ".nominals.all.chunks.txt.gz", sep="") 
   x <- fread(fn, sep=" ", header=F, data.table=F, stringsAsFactors=F)
   names(x) <- c("gene", "SNP", "DTSS", "pval_vQTL", "beta_vQTL")
    
   x2 <- x%>%mutate(gene_SNP=paste(gene, SNP, sep="_"))%>%
       dplyr::filter(gene_SNP%in%res2$gene_SNP)%>%
       mutate(zscore_vQTL=abs(qnorm(pval_vQTL*0.5))*sign(beta_vQTL),
              SE_vQTL=abs(beta_vQTL/zscore_vQTL))
   x2 <- x2[, c("gene_SNP", "beta_vQTL", "SE_vQTL", "zscore_vQTL", "pval_vQTL")] 

   ##
   res2 <- res2%>%left_join(x2, by="gene_SNP")
   res2
})

opfn <- "./5_summary.outs/dap-g_pct_0.02_union/vQTL_annoted_comb.xlsx"
write.xlsx(resNew, file=opfn, overwrite=T)




######################################
### obtain annotation information  ###
######################################

###
###

## fn <- "/nfs/rprdata/julong/sc-atac/analyses.2021-02-05/3_motif/2_motif.activities.outs/7.3_NKcell_union.motif.txt"
## resp_motif <- unique(read.table(fn, header=T)$gene)


prefix <- "/nfs/rprdata/julong/sc-atac/genetic.analysis_torus_2021-10-14/SNPannotation/"
fn <- paste(prefix, "motifList2020.txt", sep="")
motifList <- read.table(fn)$V1
#### motif 
snpmotif <- lapply(motifList, function(motif){
###
## cat(motif, "\n") 
    fn <- paste(prefix, "annot_jaspar2020/allsnp_",  motif, ".bed.gz", sep="")   
   xx <- try(fread(fn, header=F, data.table=F, stringsAsFactors=F), silent=T)
   if ( class(xx)!="try-error"){ 
      xx$motifs <- motif
   }else{
      xx <- NULL
   }   
   xx
})

snpmotif <- do.call(rbind, snpmotif[!is.na(snpmotif)])
snpmotif <- snpmotif%>%mutate(chr_pos=paste(V1, V2, sep="_"))

opfn <- gzfile(paste0(prefix, "annot_jaspar2020/zzz_allmotifs.bed.gz"))
write.table(snpmotif, opfn, row.names=F, quote=F)



###
###
fn <- "../../genetic.analysis_torus_2021-10-14/SNPannotation/annot_jaspar2020/zzz_allmotifs.bed.gz"
snpmotif <- fread(fn, header=T, data.table=F)


###
###
fn <- "./5_summary.outs/dap-g_pct_0.02_union/vQTL_annoted_summary.xlsx"
res2 <- read.xlsx(fn)
res3 <- res2%>%filter(peaking_d==3)

##
fn <- "/nfs/rprdata/julong/sc-atac/analyses.2021-02-05/3_motif/2_motif.activities.outs/3_motif.diff.results.rds"
resDiff <- read_rds(fn)




MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
resp_motif <- unique(read.table(fn, header=T)$gene)
res_anno <- NULL
for (i in 1:nrow(res3)){
   ###
   oneMCl <- res3$MCls[i]
   ii <- which(MCls==oneMCl)
   prefix <- "/nfs/rprdata/julong/sc-atac/analyses.2021-02-05/3_motif/2_motif.activities.outs/" 
   fn <- paste0(prefix, "7.", ii, "_", oneMCl, "_union.motif.txt", sep="")
   resp_motif <- read.table(fn, header=T)$gene
    
   ### 
   pos <- res3$chr_pos[i]    
   motif0 <- snpmotif%>%
      filter(chr_pos==pos)%>%pull(motifs)
   motif2 <- intersect(resp_motif, motif0)
   ###
    
   x <- resDiff%>%filter(MCls==oneMCl, gene%in%motif2) 
   x2 <- x%>%dplyr::select(MCls, contrast, motif, beta, pval, qval)

   res4 <- res3[i, c("conditions", "gene", "SYMBOL", "SNPs", "pip")]  
   x2 <- cbind(x2, res4)
   res_anno <- rbind(res_anno, x2)
}    

opfn <- "./5_summary.outs/dap-g_pct_0.02_union/vQTL_resp_motif.xlsx"
write.xlsx(res_anno, opfn, overwrite=T)



###
### genotypes
## obtain vcf files
fn <- "./5_summary.outs/dap-g_pct_0.02_union/vQTL_annoted_summary.xlsx"
res <- read.xlsx(fn)%>%filter(peaking_d==3)
chr_pos <- res$chr_pos
map_df <- str_split(chr_pos, "_", simplify=T)
##
opfn <- "./5_summary.outs/tmp_vcf/snpinfor.txt"
write.table(map_df, opfn, row.names=F, col.names=F, quote=F, sep="\t")

###
fn <- "./5_summary.outs/tmp_vcf/tmp_snp.txt"
vcf <- read.table(fn) 
gen <- vcf[,-c(1:3)]
gen <- round(gen)

DF <- NULL
for (i in 1:nrow(gen)){
   ##
   x <- as.factor(unlist(gen[i,]))
   ##
   df2 <- data.frame(snp_id=vcf$V3[i], "ref"=sum(x==0), "hete"=sum(x==1), "alt"=sum(x==2)) 
   DF <- rbind(DF, df2)
}    
    

##
res2 <- res%>%dplyr::select(conditions, gene, SYMBOL, snp_id=SNPs, pip)
res2 <- res2%>%left_join(DF, by="snp_id")

opfn <- "5_summary.outs/dap-g_pct_0.02_union/vQTL_anno3_infor.xlsx"
write.xlsx(res2, opfn, overwrite=T)





############################
### visulize effect size ###
############################


fn <- "5_summary.outs/dap-g_pct_0.02_union/vQTL_anno3_infor.xlsx"
res <- read.xlsx(fn)

###
plotDF <- NULL
for (i in 1:nrow(res)){
   ##
   condition <- res$conditions[i]
   geneID <- res$gene[i]
   snpID <- res$snp_id[i]
    
   cat(condition, geneID, "\n") 
    
   fn <- paste("./eQTL_results/", condition, ".nominals.all.chunks.txt.gz", sep="") 
   x <- fread(fn, sep=" ", header=F, data.table=F, stringsAsFactors=F)
   names(x) <- c("gene", "SNP", "DTSS", "pval", "beta")
    
   x2 <- x%>%dplyr::filter(gene==geneID, SNP==snpID)%>%
       mutate(zscore=abs(qnorm(pval*0.5))*sign(beta),
              SE=abs(beta/zscore))
   x2 <- x2[, c("beta", "SE", "zscore", "pval")] 

   ##
   res2 <- res[i,]
   res2 <- cbind(res2, x2)

   plotDF <- rbind(plotDF, res2)
}

plotDF <- plotDF%>%mutate(CI_upper=beta+1.96*SE, CI_lower=beta-1.96*SE)

opfn <- "./5_summary.outs/dap-g_pct_0.02_union/vQTL_anno3_infor.xlsx"
write.xlsx(plotDF, opfn, overwrite=T)



###
### snp id
## library(BSgenome.Hsapiens.UCSC.hg19)
## library(GenomicRanges)
## genome <- BSgenome.Hsapiens.UCSC.hg19
## seqlevelsStyle(genome) <- "NCBI"



###
### plots

fn <- "./5_summary.outs/dap-g_pct_0.02_union/vQTL_anno3_infor.xlsx"
plotDF <- read.xlsx(plotDF)
 
plotDF2 <- plotDF%>%dplyr::select(conditions, SYMBOL, snp_id, beta, CI_upper, CI_lower)%>%
    mutate(rs=ifelse(grepl("rs", snp_id), gsub(".*rs", "rs", snp_id), gsub(":G:A", "", snp_id)),
           gene_SNP=paste(SYMBOL, rs, sep="_"),
           condition2=gsub("-", "+", gsub("-EtOH", "", conditions)),
           ylab_axis=paste(gene_SNP, " (", condition2, ")", sep=""), 
           ylab_axis=fct_reorder(ylab_axis, condition2),
           treat=gsub(".*_", "", condition2))

col2 <-   c("CTRL"="#828282", "LPS"="#fb9a99", "LPS+DEX"="#e31a1c", "PHA"="#a6cee3", "PHA+DEX"="#1f78b4")
 
p <- ggplot(plotDF2, aes(x=beta, y=ylab_axis, color=treat))+
    geom_errorbarh(aes(xmax=CI_upper, xmin=CI_lower), size=0.5, height=0.2)+
    geom_point(shape=19, size=0.5)+
    geom_vline(aes(xintercept=0), size=0.25, linetype="dashed")+
    scale_x_continuous("Genetic effects on variability")+
    scale_color_manual(values=col2)+
    ## scale_y_discrete(sec.axis=dup_axis(~condition2))+
    theme_bw()+
    theme(legend.position="none",
          axis.title.x=element_text(size=10),
          axis.title.y=element_blank(),
          axis.text.x=element_text(size=10),
          axis.text.y=element_text(size=8))
          ##axis.title.x=element_text(size=10),
          ## axis.text.x=element_text(size=9, angle=-90, hjust=0, vjust=0.5),
          # axis.text.x=element_text(size=9),
          ##axis.text.y=element_text(size=9))

###
### 
figfn <- "./5_summary.outs/dap-g_pct_0.02_union/Fig1_vQTL_forest.png"
png(figfn, width=480, height=350, res=120)
p
dev.off()
