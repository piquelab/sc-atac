##
library(colorspace) ##, lib.loc="/wsu/el7/groups/piquelab/R/4.1.0/lib64/R/library")
library(tidyverse) ###, lib.loc="/wsu/el7/groups/piquelab/R/4.1.0/lib64/R/library")
library(data.table)
library(annotables)
library(INTACT, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
library(openxlsx)

rm(list=ls())


##    
outdir2 <- "./4_INTACT.outs/"
if ( !file.exists(outdir2)) dir.create(outdir2, recursive=T)


### INTACT for TWAS and colocalization from GTEx eQTL


#####################################################
### run INTACT to combine colocalzaition and TWAS ###
#####################################################

autosome <- as.character(1:22) 
grch38_unq <- grch38%>%
    dplyr::filter(chr%in%autosome, grepl("protein", biotype))%>%
    distinct(ensgene, chr, .keep_all=T)%>%dplyr::select(gene=ensgene, chr, biotype, symbol)

###
### filtering out protein genes


traits <- sort(read.table("traits.txt")$V1)[3]

for ( trait in traits){

cat(trait,",")
    
##    
fn <- paste("../3_enloc_gtex/enloc_output/", trait, ".enloc.gene.out", sep="")
enloc <- read.table(fn, header=T)

    
###
### smr
fn <- paste("./twas_smr.outs/", trait, "_gtex_topPIP_union_twas.txt.gz", sep="")
twas <- fread(fn, header=T, data.table=F)
names(twas)[1] <- "Gene"
 
## combine enloc and z-score of twas
df_comb <- enloc%>%inner_join(twas, by="Gene")%>%drop_na(zscore_gwas)
    
    
    
## intact analysis
res_intact <- intact(GLCP=df_comb$GLCP, z_vec=df_comb$zscore_gwas)
df_comb$PCG <- res_intact
### FDR 
df_comb <- df_comb%>%mutate(LFDR=1-PCG)%>%arrange(LFDR)
x <- df_comb$LFDR
FDR <- cumsum(x)/1:length(x)
df_comb$FDR <- FDR

df_comb <- df_comb%>%filter(Gene%in%grch38_unq$gene)

    
### 
opfn <- paste(outdir2, trait,"_gtex_topPIP_union_intact.txt", sep="")
write.table(df_comb, opfn, row.names=F, quote=F, col.names=T)

cat("INTACT genes:", sum(df_comb$FDR<0.1), "\n")    
    
}
 
### END loop for traits



################
### minimum ####
################

trait <- sort(read.table("traits.txt")$V1)[3]
    
##    
fn <- paste("../3_enloc_gtex/enloc_dtss_output/", trait, ".enloc.gene.out", sep="")

enloc <- read.table(fn, header=T)

    
###
### smr
fn <- paste("./twas_smr.outs/", trait, "_gtex_minP_twas.txt.gz", sep="")
twas <- fread(fn, header=T, data.table=F)
names(twas)[1] <- "Gene"
 
## combine enloc and z-score of twas
df_comb <- enloc%>%inner_join(twas, by="Gene")%>%drop_na(zscore_gwas)
    
    
    
## intact analysis
res_intact <- intact(GLCP=df_comb$GLCP, z_vec=df_comb$zscore_gwas)
df_comb$PCG <- res_intact
### FDR 
df_comb <- df_comb%>%mutate(LFDR=1-PCG)%>%arrange(LFDR)
x <- df_comb$LFDR
FDR <- cumsum(x)/1:length(x)
df_comb$FDR <- FDR

df_comb <- df_comb%>%filter(Gene%in%grch38_unq$gene)

    
### 
opfn <- paste(outdir2, trait,"_gtex_minP_intact.txt", sep="")
write.table(df_comb, opfn, row.names=F, quote=F, col.names=T)

cat("INTACT genes:", sum(df_comb$FDR<0.1), "\n")    
    








##########################################
### Add genetic variants annotation 
############################################



autosome <- as.character(1:22) 
geneDF <- grch38%>%
    dplyr::filter(chr%in%autosome, grepl("protein", biotype))%>%
    distinct(ensgene, chr, .keep_all=T)%>%dplyr::select(Gene=ensgene, symbol)


fn <- "SCAIP_final_bed.gz"
bed <- fread(fn, header=T, data.table=F)
bed <- bed[,c("genetic_variant", "chr_pos_grch38")]
bed <- bed%>%distinct(chr_pos_grch38, .keep_all=T)
b38 <- bed$chr_pos_grch38
names(b38) <- bed$genetic_variant

### INTACT files
trait <- sort(read.table("traits.txt")$V1)[3]
fn <- paste(outdir2, trait, "_gtex_topPIP_union_intact.txt", sep="")
res <- read.table(fn, header=T)


## torus file 
fn <- "/nfs/rprdata/julong/sc-atac/genetic.analysis_torus_2021-10-14/SNPannotation/4_SNPAnnot.outs/pct_0.02_combineNew/Union_torus.annot.gz"
anno <- fread(fn, header=T, data.table=F)
names(anno)[1] <- "genetic_variant"
anno$chr_pos_grch38 <- b38[anno$genetic_variant]
anno <- anno%>%drop_na(chr_pos_grch38)

peaking <- anno$peaking_d
names(peaking) <- anno$chr_pos_grch38

### twas file
fn <- paste("./3_twas_summary.outs/", trait, "_gtex_topPIP_union_twas.txt.gz", sep="")
twas <- fread(fn, header=T, data.table=F)
twas2 <- twas%>%dplyr::select(Gene=gene, FDR_twas=FDR)


res2 <- res%>%left_join(geneDF, by="Gene")%>%left_join(twas2, by="Gene")

res2 <- res2%>%mutate(peaking_d=ifelse(chr_pos_grch38%in%anno$chr_pos_grch38, peaking[chr_pos_grch38], 0),
                      comb=paste(Gene, genetic_variant, sep="_"))




fn <- "/nfs/rprdata/julong/sc-atac/twas_analysis_2022-10-16/Whole_Blood_GTEx_v8_results/Whole_Blood.fastqtl.gz"
fast <- fread(fn, header=F, data.table=F)
names(fast) <- c("gene", "genetic_variant", "DTSS", "pval", "beta", "SE")
fast <- fast%>%mutate(gene=gsub("\\..*", "", gene))

pval <- fast$pval
names(pval) <- paste(fast$gene, fast$genetic_variant, sep="_")

### Add pval_eqtl 
res2$pval_eqtl <- pval[as.character(res2$comb)]

## reorder column name
res2 <- res2%>%
    dplyr::select(Gene, symbol, genetic_variant, chromosome=chr, pos, chr_pos_grch38, PIP, pval_eqtl,
                  peaking_d, pval_twas=pval_gwas, zscore_twas=zscore_gwas, FDR_twas, GLCP, PCG, FDR_intact=FDR)

res2 <- res2%>%arrange(FDR_intact)

opfn <- paste(outdir2, trait, "_gtex_topPIP_union_combinfor.txt", sep="")
write.table(res2, file=opfn, quote=F, row.names=F, sep="\t")



###
### add chr_pos_grch37

fn <- "SCAIP_final_bed.gz"
bed <- fread(fn, header=T, data.table=F)
bed <- bed[,c("chr_pos", "chr_pos_grch38")]
bed <- bed%>%distinct(chr_pos_grch38, .keep_all=T)
b37 <- bed$chr_pos
names(b37) <- bed$chr_pos_grch38


trait <- sort(read.table("traits.txt")$V1)[3]
fn <- paste(outdir2, trait, "_gtex_topPIP_union_combinfor.txt", sep="")
res <- fread(fn, header=T, data.table=F, sep="\t")

res <- res%>%mutate(chr_pos_grch37=b37[chr_pos_grch38])

res2 <- res%>%
    dplyr::select(Gene, symbol, genetic_variant, chromosome, pos, chr_pos_grch38, chr_pos_grch37,
                  PIP, pval_eqtl,
                  peaking_d, pval_twas, zscore_twas, FDR_twas, GLCP, PCG, FDR_intact)
###
###
opfn <- paste(outdir2, trait, "_gtex_topPIP_union_combinfor.txt", sep="")
write.table(res2, file=opfn, quote=F, row.names=F, sep="\t")




##################################################
### compare annotation and without annotation 
##################################################

trait <- sort(read.table("traits.txt")$V1)[3]


###
### aloft, topPIP with annotation

fn <- "./5_pub.outs/2_supp_tables/TableS5_1_asthma-risk-genes_ALOFT.txt.gz"
res <- read.table(fn, header=T, sep="\t")
sig <- res%>%filter(FDR_intact<0.1)%>%pull(Gene)%>%unique()


###
### gtex, topPIP, annotation
fn <- "./5_pub.outs/2_supp_tables/TableS6_1_asthma-risk-genes_gtex.txt.gz"
res_gtex <- read.table(fn, header=T, sep="\t")
sig_gtex <- res_gtex%>%filter(FDR_intact<0.1)%>%pull(Gene)%>%unique()

sig <- union(sig, sig_gtex)


###
### minP w/o annotation

fn <- "./4_INTACT.outs/Asthma_Bothsex_inv_var_meta_GBMI_052021_nbbkgt1_aloft_minP_intact.txt"
res2 <- read.table(fn, header=T)
sig2 <- res2%>%filter(FDR<0.1)%>%pull(Gene)%>%unique()

fn <- "./4_INTACT.outs/Asthma_Bothsex_inv_var_meta_GBMI_052021_nbbkgt1_gtex_minP_intact.txt"
res2_gtex <- read.table(fn, header=T)
sig2_gtex <- res2_gtex%>%filter(FDR<0.1)%>%pull(Gene)%>%unique()


sig2 <- union(sig2, sig2_gtex)


shared <- intersect(sig, sig2)
unq <- setdiff(sig, sig2)



















