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



#####################################################
### run INTACT to combine colocalzaition and TWAS ###
#####################################################

autosome <- as.character(1:22) 
grch38_unq <- grch38%>%
    dplyr::filter(chr%in%autosome, grepl("protein", biotype))%>%
    distinct(ensgene, chr, .keep_all=T)%>%dplyr::select(gene=ensgene, chr, biotype, symbol)

###
### filtering out protein genes


traits <- sort(read.table("traits.txt")$V1)[3:4]

for ( trait in traits){

cat(trait,",")
    
##    
fn <- paste("../3_enloc_aloft/enloc_output/", trait, ".enloc.gene.out", sep="")
enloc <- read.table(fn, header=T)

    
###
### smr
fn <- paste("./twas_smr.outs/", trait, "_aloft_topPIP_union_twas.txt.gz", sep="")
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
opfn <- paste(outdir2, trait,"_aloft_topPIP_union_intact.txt", sep="")
write.table(df_comb, opfn, row.names=F, quote=F, col.names=T)

cat("INTACT genes:", sum(df_comb$FDR<0.1), "\n")    
    
}
 
### END loop for traits


#######################################################################
### twas from minimum pvalue and colocalization also use annotation 
########################################################################

trait <- sort(read.table("traits.txt")$V1)[3]
    
##    
fn <- paste("../3_enloc_aloft/enloc_dtss_output/", trait, ".enloc.gene.out", sep="")
enloc <- read.table(fn, header=T)

    
###
### smr
fn <- paste("./twas_smr.outs/", trait, "_aloft_minP_twas.txt.gz", sep="")
twas <- fread(fn, header=T, data.table=F)
names(twas)[1] <- "Gene"

twas <- twas%>%dplyr::select(Gene, genetic_variant, pval_eqtl=pval, chr_pos_grch38, chr, pval_gwas, zscore_gwas)
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
opfn <- paste(outdir2, trait,"_aloft_minP_intact.txt", sep="")
write.table(df_comb, opfn, row.names=F, quote=F, col.names=T)



##################################################
### compare annotation and without annotation 
##################################################

trait <- sort(read.table("traits.txt")$V1)[3]

###
### topPIP, annotation
fn <- paste(outdir2, trait, "_aloft_topPIP_union_intact.txt", sep="")
res <- read.table(fn, header=T)
sig <- res%>%filter(FDR<0.1)%>%pull(Gene)%>%unique()

###
### minP, w/o annotation
fn <- paste(outdir2, trait, "_aloft_minP_intact.txt", sep="")
res2 <- read.table(fn, header=T)
sig2 <- res2%>%filter(FDR<0.1)%>%pull(Gene)%>%unique()

olap <- intersect(sig, sig2)







##########################################
### Add genetic variants annotation 
############################################

autosome <- as.character(1:22) 
grch38_unq <- grch38%>%
    dplyr::filter(chr%in%autosome, grepl("protein", biotype))%>%
    distinct(ensgene, chr, .keep_all=T)%>%dplyr::select(gene=ensgene, chr, biotype, symbol)



###
geneDF <- grch38_unq%>%dplyr::select(Gene=gene, symbol)

### INTACT files
trait <- sort(read.table("traits.txt")$V1)[3]
fn <- paste(outdir2, trait, "_aloft_topPIP_union_intact.txt", sep="")
res <- read.table(fn, header=T)

## torus file
fn <- "/nfs/rprdata/julong/sc-atac/genetic.analysis_torus_2021-10-14/SNPannotation/4_SNPAnnot.outs/pct_0.02_combineNew/Union_torus.annot.gz"
anno <- fread(fn, header=T, data.table=F)
names(anno)[1] <- "genetic_variant"

### twas file
fn <- paste("./3_twas_summary.outs/", trait, "_aloft_topPIP_union_twas.txt.gz", sep="")
twas <- fread(fn, header=T, data.table=F)
twas2 <- twas%>%dplyr::select(Gene=gene, FDR_twas=FDR)

res2 <- res%>%left_join(anno, by="genetic_variant")%>%left_join(geneDF, by="Gene")%>%left_join(twas2, by="Gene")

opfn <- paste(outdir2, trait, "_aloft_topPIP_union_combinfor.txt", sep="")
write.table(res2, file=opfn, quote=F, row.names=F, sep="\t")



##############################
### combine eqtl results
##############################

###
### summary data
trait <- sort(read.table("traits.txt")$V1)[3]
fn <- paste(outdir2, trait, "_aloft_topPIP_union_combinfor.txt", sep="")
res <- read.table(fn, header=T, sep="\t")
res <- res%>%mutate(comb=paste(Gene, genetic_variant, sep="_"))


###
### fast eqtl results
fn <- "/nfs/rprdata/julong/sc-atac/twas_analysis_2022-10-16/ALOFT_results/PC1-18.nominals.eQTL.txt.gz"
fast <- fread(fn, header=F, data.table=F)
names(fast) <- c("gene", "genetic_variant", "DTSS", "pval", "beta")

pval <- fast$pval
names(pval) <- paste(fast$gene, fast$genetic_variant, sep="_")

### Add pval_eqtl 
res$pval_eqtl <- pval[as.character(res$comb)]


### output
opfn <- paste(outdir2, trait, "_aloft_topPIP_union_combinfor.txt", sep="")
write.table(res, file=opfn, quote=F, row.names=F, sep="\t")



## nsig <- res2%>%filter(FDR<0.1)%>%pull(symbol)%>%length()
## nsig2 <- res2%>%filter(FDR<0.1, peaking_d==3)%>%pull(symbol)%>%length()

### bed file
## fn <- "SCAIP_final_bed.gz"
## bed <- fread(fn, header=T, data.table=F)





###############################
### check the script 
###############################


####################################
### keep the same to old 
##################################


## autosome <- as.character(1:22) 
## grch38_unq <- grch38%>%
##     dplyr::filter(chr%in%autosome, grepl("protein", biotype))%>%
##     distinct(ensgene, chr, .keep_all=T)%>%dplyr::select(gene=ensgene, chr, biotype, symbol)



## trait <- sort(read.table("traits.txt")$V1)[4]
    
## ##    
## fn <- paste("../3_enloc_aloft/enloc_output/", trait, "_old.enloc.gene.out", sep="")
## enloc <- read.table(fn, header=T)

    
## ###
## ### smr
## fn <- paste("./twas_smr.outs/", trait, "_aloft_topPIP_union_twas.txt.gz", sep="")
## twas <- fread(fn, header=T, data.table=F)
## names(twas)[1] <- "Gene"
 
## ## combine enloc and z-score of twas
## df_comb <- enloc%>%inner_join(twas, by="Gene")%>%drop_na(zscore_gwas)%>%
##     filter(Gene%in%grch38_unq$gene)
    
    
## ## intact analysis
## res_intact <- intact(GLCP=df_comb$GLCP, z_vec=df_comb$zscore_gwas)
## df_comb$PCG <- res_intact
## ### FDR 
## df_comb <- df_comb%>%mutate(LFDR=1-PCG)%>%arrange(LFDR)
## x <- df_comb$LFDR
## FDR <- cumsum(x)/1:length(x)
## df_comb$FDR <- FDR

## ## df_comb <- df_comb

    
## ###  
## opfn <- paste(outdir2, trait,"_old_aloft_topPIP_union_intact.txt", sep="")
## write.table(df_comb, opfn, row.names=F, quote=F, col.names=T)

