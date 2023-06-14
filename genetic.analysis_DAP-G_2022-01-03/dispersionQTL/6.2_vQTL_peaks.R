###
library(tidyverse)
library(Seurat)
library(Signac)
###
library(ChIPseeker, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
library(cowplot)
library(openxlsx)
library(ggrepel)

###
rm(list=ls())


outdir <- "./6_vQTL_summ.outs/vQTLs_peak/"
if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)

###
### evaluate the effects of genetic variants on peaks 


#########################
### vQTLs information ###
#########################

fn <- "./5_summary.outs/dap-g_pct_0.02_union/vQTL_annoted_comb.xlsx"
res <- read.xlsx(fn)
res2 <- res%>%filter(peaking_d>0)


#############################
### obtain genotype files ###
#############################

chr_pos <- unique(res$chr_pos)
df <- str_split(chr_pos, "_", simplify=T)
##
opfn <- "./6_vQTL_summ.outs/vcf/snpinfor.txt" 
write.table(df, opfn, row.names=F, col.names=F, quote=F, sep="\t")

genfn <- "./6_vQTL_summ.outs/vcf/vQTLs_dosages.txt"
vcf <- read.table(genfn, header=F, fill=T)
##
sampleID <- read.table("./6_vQTL_summ.outs/vcf/sample6.id")$V1
gen <- vcf[,-c(1:3)]
rownames(gen) <- vcf$V3
colnames(gen) <- sampleID



#############################
### get pseudo-bulk reads ###
#############################

fn <- "../../analyses.2021-02-05/2_Differential/1.2_DiffPeak.outs/1_YtX.sel.rds"
mat <- as.data.frame(read_rds(fn))
### normalized data
depth <- 1/colSums(mat)

### sample information
cvt <- str_split(colnames(mat), "_", simplify=T)
meta <- data.frame(rn=colnames(mat), conditions=paste(cvt[,1], cvt[,2], sep="_"), sampleID=cvt[,3])


## peak information
DFpeak <- as.data.frame(str_split(rownames(mat), "-", simplify=T))
names(DFpeak) <- c("chr", "start", "end") 
DFpeak$peak <- rownames(mat)



###
### individual level data for plots 
plotDF <- NULL
for (i in 1:nrow(res2)){
###    
ii <- gsub("-EtOH", "", res2$conditions[i])    
    
snp_i <- res2$SNPs[i]
chr_i <- res2$chr[i]
pos_i <- res2$pos[i]
    
##    
meta2 <- meta%>%filter(conditions==ii)    

    
###
### genotype
gen_i <- unlist(gen[snp_i, meta2$sampleID])    
meta2$dosage <- gen_i
    
###    
### peak reads    
peakSel <- DFpeak%>%filter(chr==chr_i, pos_i>=as.numeric(start), pos_i<=as.numeric(end))%>%pull(peak)  
ysum <- colSums(mat[peakSel, meta2$rn])/depth[meta2$rn]
ysum2 <- log2(ysum*1e+06+1)
meta2$y <- ysum2

meta2 <- meta2%>%
   mutate(gene=res2$gene[i], symbol=res2$SYMBOL[i], SNP=snp_i, snpinfor=res2$snpinfor[i],
          peaking_d=res2$peaking_d[i]) 

###                        
plotDF <- rbind(plotDF, meta2)
cat(i, length(peakSel), "\n")       
}

### output
opfn <- paste(outdir, "1_vQTLs_individual.rds", sep="")
write_rds(plotDF, file=opfn)


######################################
### Regression genotype with peaks ###
######################################

resNew <- NULL
for (i in 1:nrow(res2)){
###    
ii <- gsub("-EtOH", "", res2$conditions[i])    
    
snp_i <- res2$SNPs[i]
chr_i <- res2$chr[i]
pos_i <- res2$pos[i]

##    
meta2 <- meta%>%filter(conditions==ii)    
    
###
### genotype
gen_i <- unlist(gen[snp_i, meta2$sampleID])    
meta2$dosage <- gen_i
    
###    
### peak reads    
peakSel <- DFpeak%>%filter(chr==chr_i, pos_i>=as.numeric(start), pos_i<=as.numeric(end))%>%pull(peak)  
ysum <- colSums(mat[peakSel, meta2$rn])/depth[meta2$rn]
ysum2 <- log2(ysum*1e+06+1)
meta2$y <- ysum2
nind <- nrow(meta2)

###    
### effects of genetic variants on peaks    
lm0 <- lm(y~dosage, data=meta2)    
beta_peak <- coef(lm0)["dosage"]
SE_peak <- sqrt(vcov(lm0)["dosage", "dosage"])    
zscore_peak <- beta_peak/SE_peak
pval_peak <- pt(abs(zscore_peak), df=nind-2, lower.tail=F)*2    
    
df <- data.frame(beta_peak, SE_peak, zscore_peak, pval_peak)

    
tmp <- cbind(res2[i,], df)
resNew <- rbind(resNew, tmp)
    
}

opfn <- paste(outdir, "2_vQTLs_comb_infor.xlsx", sep="")
write.xlsx(resNew, opfn, overwrite=T)

### short for results
opfn <- paste(outdir, "2.2_vQTLs_comb.xlsx", sep="")
x <- resNew%>%dplyr::select(conditions, gene, SNPs, SYMBOL,pip,  peaking_d,
   beta_vQTL, SE_vQTL, zscore_vQTL, pval_vQTL, beta_peak, SE_peak, zscore_peak, pval_peak)
write.xlsx(x, file=opfn, overwrite=T)



####
### scatter plots
fn <- paste(outdir, "2.2_vQTLs_comb.xlsx", sep="")
plotDF <- read.xlsx(fn)
plotDF2 <- plotDF[-c(14,15),]

###
x <- plotDF2$zscore_peak
y <- plotDF2$zscore_vQTL
corr0 <- cor.test(x, y, method="pearson")
rr <- round(corr0$estimate, 3)
lab1 <- bquote(~"PCC"==.(rr)~"NS")
##
corr1 <- cor.test(x, y, method="spearman")
rr <- round(corr1$estimate, 3)
lab2 <- bquote(~rho==.(rr)~"*")

p <- ggplot(plotDF2, aes(x=zscore_peak, y=zscore_vQTL))+
   geom_point(shape=24, size=2.5)+
   geom_text_repel(aes(x=zscore_peak, y=zscore_vQTL, label=SYMBOL), color="grey30",
                   size=2.5, fontface="italic", max.overlaps=20)+
   annotate("text", x=-2, y=-1, label=lab1, size=3, color="grey30", hjust=0)+
   annotate("text", x=-2, y=-2, label=lab2, size=3, color="grey30", hjust=0)+ 
   xlab(bquote(~italic(Z)~"-score of genetic variants on peak"))+ 
   ylab(bquote(~italic(Z)~"-score of genetic variants on variability"))+
   geom_hline(yintercept=0, linetype="dashed", color="grey30")+
   geom_vline(xintercept=0, linetype="dashed", color="grey30")+ 
   theme_bw()+
   theme(axis.text=element_text(size=10),
         axis.title=element_text(size=10))

###
figfn <- paste(outdir, "Figure1_scatter_zscore.png", sep="")
png(figfn, width=420, height=420, res=120)
p
dev.off()


       

####################
### violin plots ###
####################

fn <- paste(outdir, "1_vQTLs_individual.rds", sep="")
plotDF <- read_rds(fn)


fn <- paste(outdir, "2_vQTLs_comb_infor.xlsx", sep="")
vQTL <- read.xlsx(fn)%>%dplyr::filter(peaking_d==3|SYMBOL%in%c("TNFAIP2", "MYCBP2"))
vQTL <- vQTL%>%arrange(SYMBOL)

###
### plots

figs_ls <- lapply(1:nrow(summ), function(i){
###    
condition <- gsub("-EtOH", "", vQTL$conditions[i])
geneID <- vQTL$gene[i]
symbol_i <- vQTL$SYMBOL[i]
snpID <- vQTL$SNPs[i]

cat(symbol_i, condition, snpID, "\n")

### plot data    
plotDF2 <- DF2%>%dplyr::filter(conditions==condition, symbol==symbol_i)%>%
    mutate(genotype=as.character(round(dosage, 0))) 

    
### plot setting
snpID2 <- gsub(";.*", "", snpID)
snp_ls <- unlist(str_split(snpID2, ":"))
ref <- snp_ls[3]
alt <- snp_ls[4]
gen_lab <- c("0"=paste0(ref, ref), "1"=paste0(ref, alt), "2"=paste0(alt, alt))
gen_alpha <- c("0"=0.2, "1"=0.6, "2"=1)


### plots
condition2 <- gsub("-", "+", condition)    
p2 <- ggplot(plotDF2, aes(x=factor(genotype), y=y))+
   geom_violin(aes(alpha=factor(genotype)), fill="grey40", width=0.8, lwd=0.2)+
   ##geom_boxplot(width=0.2, color="grey", outlier.shape=NA)+
   geom_jitter(color="red", width=0.2, size=0.6)+ 
   scale_alpha_manual(values=gen_alpha)+
   scale_x_discrete(labels=gen_lab)+
   ylab(bquote("Normalized peaks"~"("~log[2]~"CPM"~")"))+
   ##geom_smooth(aes(x=dosage, y=y), method="lm", se=F, color="blue")+ 
   ggtitle(bquote(~italic(.(symbol_i))~"in"~.(condition2)))+ 
   theme_bw()+
   theme(legend.position="none",
         axis.text=element_text(size=10),
         axis.title.x=element_blank(),
         axis.title.y=element_text(size=10),
         plot.title=element_text(size=10, hjust=0.5))
## figfn <- paste(outdir, "Figure2_", i, "_", geneID, "_", symbol_i, "_peak.violin.png", sep="")
## png(figfn, width=320, height=320, res=120)
## print(p2)
## dev.off()
p2    
})

### output
figfn <- paste(outdir, "Figure2_peak_violin.png", sep="")
png(figfn, width=1100, height=550, res=120)
plot_grid(plotlist=figs_ls, ncol=4)
dev.off()
