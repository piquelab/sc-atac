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


outdir <- "./5_pub.outs/3_example_plots/"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)



#################################
### Obtain genotype vcf files ###
#################################

 
outdir2 <- paste(outdir, "vcf/", sep="")
if (!file.exists(outdir2)) dir.create(outdir2, showWarnings=F, recursive=T)
 
trait <- sort(read.table("traits.txt")$V1)[3]

###                            
### candidate genes 
fn <- paste("./4_INTACT.outs/", trait, "_aloft_topPIP_union_combinfor.txt", sep="")
res <- read.table(fn, header=T, sep="\t")

res2 <- res%>%filter(PIP>0.1, FDR<0.1, FDR_twas<0.1)

###
bed <- str_split(res2$genetic_variant, ":", simplify=T)
res2 <- res2%>%mutate(pos=as.integer(bed[,2]))
 
for ( i in 1:nrow(res2)){
   ##
   chr <- res2$chr[i]
   pos <- res2$pos[i]
   pos1 <- as.integer(pos-1e+06)
   pos2 <- as.integer(pos+1e+06)

   ens <- res2$Gene[i]
   symbol <- res2$symbol[i]
    
   ### vcf files
   genfn <- paste(outdir2, ens, "_", symbol, "_", chr, "_dosages.txt", sep="")    
   system(paste("bcftools query -r ", chr, ":", pos1, "-", pos2,
                " -f '%CHROM\\t%POS\\t%ID[\\t%DS]\\n' ../../ALOFT.vcf.gz > ", genfn, sep=""))
   ##
   cat(symbol, "\n")
}    



##############################
### locuszoom plots 
##############################



#################
### Read data ###
#################


###
### grch37 corresponding to grch38
bed <- fread("SCAIP_final_bed.gz", header=T, data.table=F)
bed2 <- bed%>%distinct(chr_pos_grch38, .keep_all=T)

## names(pos38) <- bed2$genetic_variant

pos37 <- bed2$chr_pos
names(pos37) <- bed2$chr_pos_grch38



###
### candidate genes
trait <- sort(read.table("traits.txt")$V1)[3]
fn <- paste("./4_INTACT.outs/", trait, "_aloft_topPIP_union_combinfor.txt", sep="")
res <- read.table(fn, header=T, sep="\t")

x <- str_split(res$genetic_variant, ":", simplify=T)
res <- res%>%mutate(chr_pos_grch37=paste(x[,1], x[,2], sep="_"))




###
### gwas 
fn <- paste("./gwas_impute/", trait, "_impute.txt.gz", sep="")
gwas <- fread(fn, header=T, data.table=F)
gwas <- gwas%>%mutate(chr_pos_grch37=pos37[id_b38_0])

###
### dap results 
fn <- "/nfs/rprdata/julong/sc-atac/twas_analysis_2022-10-16/ALOFT_results/aloft_allSNPs_union.txt.gz"
dap <- fread(fn, header=F, data.table=F)
names(dap) <- c("gene", "id_b37", "PIP")


###
### FastQTL results
fn <- "/nfs/rprdata/julong/sc-atac/twas_analysis_2022-10-16/ALOFT_results/PC1-18.nominals.eQTL.txt.gz"
fast <- fread(fn, data.table=F)
names(fast) <- c("gene", "genetic_variant", "DTSS", "pval", "beta")


###
### signac object
fn <- "/nfs/rprdata/julong/sc-atac/analyses.2021-02-05/1_processing/5.1_reCallPeak.outs/3_scATAC.annot.rds"
atac <- read_rds(fn)
atac <-subset(atac, MCls!="DC")
x <- atac@meta.data
x2 <- x%>%mutate(treat=gsub(".*-ATAC-|_.*", "", NEW_BARCODE))
atac <- AddMetaData(atac, metadata=x2)

###
### sample list
sampleID <- read.table("sampleID.txt")$V1

###
### differential accessible results
## fn <- paste(outdir2, "zzz2_DAR.xlsx", sep="")
## DAR <- openxlsx::read.xlsx(fn)



###
### contrast
contrast_ls <- list("LPS"=c("LPS", "CTRL"), "LPS-DEX"=c("LPS-DEX", "LPS"),
     "PHA"=c("PHA", "CTRL"), "PHA-DEX"=c("PHA-DEX", "PHA"))               

col_treat <- c("LPS"="#fb9a99", "LPS-DEX"="#e31a1c", "PHA"="#a6cee3", "PHA-DEX"="#1f78b4")


##############
#### plots ###
##############
 
## gene1 <- tmp2%>%filter(peaking_d==3)%>%pull(symbol)
## fn <- "./4.3_example.outs/ALOFT_union/genes_list3/zzz_annot.xlsx"
## DFanno <- openxlsx::read.xlsx(fn)


outdir2 <- "./5_pub.outs/3_example_plots/Examples2/"
if (!file.exists(outdir2)) dir.create(outdir2, showWarnings=F, recursive=T)

### supp
## geneSel <- c("GSDMB", "ORMDL3", "SLC22A5", "MICB", "FLOT1", "ERBB3", "WDR36", "FADS2")
## res3 <- res2%>%dplyr::filter(symbol%in%geneSel)
fn <- "./5_pub.outs/3_example_plots/TableS5_1_asthma-risk-genes_ALOFT.txt.gz"
res <- fread(fn, header=T, sep="\t", data.table=F)

## main
res2 <- res%>%dplyr::filter(FDR_intact<0.1, FDR_twas<0.1, PIP>0.1, chromosome==17)%>%arrange(FDR_twas)

## opfn <- paste(outdir, "geneList.xlsx", sep="")
## write.xlsx(res2, file=opfn)
## opfn <- paste(outdir, "geneList.txt", sep="")
## write.table(res2, file=opfn, quote=F, row.names=F)

## geneSel <- c("KIF1B", "TNFRSF14", "SIK2", "NDFIP1",  "TRAT1", "HAUS8", "MPHOSPH9")[c(1,2,4)]
## res3 <- res3%>%filter(symbol%in%geneSel)
## res3 <- res3[c(3,4,8),]
    ##filter(peaking_d==dd, PIP>0.5, FDR_intact<0.1)
res3 <- res2[1:10,]

x <- str_split(res3$genetic_variant, ":", simplify=T)
res3$chr_pos_grch37 <- paste(x[,1], x[,2], sep="_")

###############
### plots-1 ###
###############
 
##dtss <- 2.5e+04 ## default 5e+04, 10e+04
dtss1 <- 5e+04*2
dtss2 <- 5e+04*2

 

for (i in  1:10){
###    

ens <- res3$Gene[i]
symbol <- res3$symbol[i]
symbol2 <- symbol    

chr_pos <- res3$chr_pos_grch37[i]
chr_i <- as.numeric(gsub("_.*", "", chr_pos))
pos_i <- as.numeric(gsub(".*_", "",  chr_pos))
snp_i <- res3$genetic_variant[i]
    
cat(i, ens, symbol, "\n")

###if (grepl("HLA", symbol)) next

#########################
### local zoom  plots ###
#########################

    
#############################
### 0, calculate LD files ###    
#############################

###
### vcf file
genfn <- paste("./5_pub.outs/3_example_plots/vcf/", ens, "_", symbol, "_", chr_i, "_dosages.txt", sep="")
gen <- read.table(genfn, header=F)

###
gen2 <- gen%>%filter(V2>pos_i-dtss1, V2<pos_i+dtss2)
isnp <- which(gen2$V2==pos_i)

### LD data frame    
ldDF <- NULL
g0 <- as.numeric(gen2[isnp,4:ncol(gen2)])    
for (k in 1:nrow(gen2)){
   ##
   g1 <- as.numeric(gen2[k,4:ncol(gen2)])
   r2 <- cor(g0, g1)^2
   ## 
   DF <- data.frame(chr=gen2[k,1], pos=gen2[k,2], genetic_variant=gen2[k,3], r2=round(r2, digits=3))
   ldDF <- rbind(ldDF, DF) 
}    

ldDF2 <- ldDF%>%
    mutate(chr_pos_grch37=paste(chr, pos, sep="_"))%>%
    dplyr::select(chr_pos_grch37, r2)%>%
    mutate(gr=case_when(r2>=0.8~"1", r2>=0.6&r2<0.8~"2", r2>=0.4&r2<0.6~"3",
                        r2>=0.2&r2<0.4~"4", TRUE~"5"))



    
    
####################
### plot 1, gwas ###
####################

    
gwas2 <- gwas%>%filter(chr==chr_i)%>%drop_na(chr_pos_grch37)
gwas2 <- gwas2%>%
   mutate(pos0=as.integer(gsub(".*_", "", chr_pos_grch37)), pos2=pos0/1e+03, log10p=-log10(pval))%>%
    filter(pos0>(pos_i-dtss1), pos0<(pos_i+dtss2))

###
###
gwas3 <- gwas2%>%
    left_join(ldDF2, by="chr_pos_grch37")
gwas3 <- gwas3%>%mutate(gr2=ifelse(chr_pos_grch37==chr_pos, 0, gr)) ##, gr2=ifelse(r2==1&!is.na(r2), 0, gr2))
    
##
anno <- gwas3%>%filter(chr_pos_grch37==chr_pos)

    
p1 <- ggplot(gwas3, aes(x=pos2, y=log10p))+
   geom_point(aes(colour=factor(gr2), shape=factor(gr2), size=factor(gr2) ))+
   scale_size_manual("",
       values=c("0"=3.5, "1"=2,"2"=1.5, "3"=1.5, "4"=1.5, "5"=1.5, "6"=1.5), guide="none")+
   ## scale_color_manual("",
   ##     values=c("1"="#D43F3A", "2"="#EEA236", "3"="#5CB85C", "4"="#46B8DA",
   ##              "5"="#357EBD", "6"="#B8B8B8", "0"="#9632B8"), guide="none")+      
   scale_color_manual("",
       values=c("1"="#D43F3A", "2"="#EEA236", "3"="#5CB85C", "4"="#46B8DA",
                "5"="#357EBD", "6"="#B8B8B8", "0"="#9632B8"),
       breaks=c(1, 2, 3, 4, 5),
       labels=c("1"="0.8~1", "2"="0.6~0.8", "3"="0.4~0.6",
                "4"="0.2~0.4", "5"="0~0.2"),
       guide=guide_legend(override.aes=list(size=2)))+
   scale_shape_manual("",
       values=c("0"=18, "1"=19, "2"=19, "3"=19, "4"=19, "5"=19, "6"=19), guide="none")+    
   geom_text_repel(data=anno, aes(x=pos2, y=log10p, label=rs),
                   color="#9632B8", size=3.5,
                   min.segment.length=0, segment.curvature=1e-20, box.padding=0.8)+
   ## ggtitle(bquote(~italic(.(symbol)) ))+     
   ## annotate("point", x=as.numeric(pos38_i/1e+03), y=0, shape=18, size=3.5, color="blue")+ 
   xlab(bquote("position (kb)"~"(chr"~.(chr_i)~")"))+
   ylab(bquote(~"GWAS-asthma"~-log[10]~"("~italic(p)~")"))+
   theme(plot.title=element_text(hjust=0.5, size=12),
      axis.text=element_text(size=10),
      axis.title=element_text(size=10),
      legend.position="none",
      legend.key=element_blank(),
      legend.key.size=grid::unit(0.4, "cm"),
      legend.background=element_blank(),
      legend.box.background=element_blank(),
      legend.text=element_text(size=8),
      plot.margin=unit(c(5.5, 12, 5.5, 5.5), "pt"),
       panel.background=element_blank(),
       panel.border=element_rect(fill=NA,color="black"),
       panel.grid.major=element_blank(),
       panel.grid.minor=element_blank())

## figfn <- paste(outdir2, "Figure", i, ".1_", ens, "_", symbol, "_gwas.png", sep="")       
## png(figfn, width=580, height=380, res=120)
## print(p1)
## dev.off()



##############################
### plot2, FastQTL results ###
##############################

## ens <- res2$gene[i]
## symbol <- res2$symbol[i] 
## chr_pos <- res2$chr_pos_grch37[i]
## chr_i <- as.numeric(gsub("_.*", "", chr_pos))
## pos_i <- as.numeric(gsub(".*_", "",  chr_pos))

fast2 <- fast%>%filter(gene==ens)
mapDF <- str_split(fast2$genetic_variant, ":", simplify=T)
fast2$chr <- as.numeric(mapDF[,1])
fast2$pos <- as.numeric(mapDF[,2])

fast2$chr_pos_grch37 <- paste(mapDF[,1], mapDF[,2], sep="_")
fast2$rs <- gsub(".*;", "", mapDF[,4])
## fast2$chr_pos_grch38 <- pos38[fast2$genetic_variant]

    
fast3 <- fast2%>%filter(chr==chr_i, pos>(pos_i-dtss1), pos<(pos_i+dtss2))%>%
    left_join(ldDF2, by="chr_pos_grch37")
    
fast3 <- fast3%>%
    mutate(pos2=pos/1e+03, log10p=-log10(pval),
           gr2=ifelse(chr_pos_grch37==chr_pos, 0, gr))


##    
anno <- fast3%>%filter(chr_pos_grch37==chr_pos)    

p2 <- ggplot(fast3, aes(x=pos2, y=log10p))+
   geom_point(aes(colour=factor(gr2), size=factor(gr2),shape=factor(gr2)))+
   scale_size_manual("",
       values=c("0"=3.5, "1"=2,"2"=1.5, "3"=1.5, "4"=1.5, "5"=1.5, "6"=1.5), guide="none")+
   ## scale_color_manual("",
   ##     values=c("1"="#D43F3A", "2"="#EEA236", "3"="#5CB85C", "4"="#46B8DA",
   ##              "5"="#357EBD", "6"="#B8B8B8", "0"="#9632B8"), guide="none")+      
   scale_color_manual("",
       values=c("1"="#D43F3A", "2"="#EEA236", "3"="#5CB85C", "4"="#46B8DA",
                "5"="#357EBD", "6"="#B8B8B8", "0"="#9632B8"),
       breaks=c(1, 2, 3, 4, 5),
       labels=c("1"="0.8~1", "2"="0.6~0.8", "3"="0.4~0.6",
                "4"="0.2~0.4", "5"="0~0.2"),
       guide=guide_legend(override.aes=list(size=2)))+
   scale_shape_manual("",
       values=c("0"=18, "1"=19, "2"=19, "3"=19, "4"=19, "5"=19, "6"=19), guide="none")+        
   geom_text_repel(data=anno, aes(x=pos2, y=log10p, label=rs),
                   color="#9632B8",
                  min.segment.length=0, segment.curvature=1e-20, 
                   size=3.5, box.padding=0.8)+
   xlab(bquote("position (kb)"~"(chr"~.(chr_i)~")"))+
   ylab(bquote(~"FastQTL"~log[10]~"("~italic(p)~")"))+
   ## ggtitle(bquote(~italic(.(symbol)) ))+ 
   theme(plot.title=element_text(hjust=0.5, size=14),
      axis.text=element_text(size=10),
      axis.title=element_text(size=10),
      plot.margin=unit(c(5.5, 12, 5.5, 5.5), "pt"),
      legend.position="none",
       panel.background=element_blank(),
       panel.border=element_rect(fill=NA,color="black"),
       panel.grid.major=element_blank(),
       panel.grid.minor=element_blank())

## figfn <- paste(outdir2, "Figure", i, ".3_", ens, "_", symbol, "_eQTL.png", sep="")       
## png(figfn, width=450, height=380, res=120)
## print(p3)
## dev.off()


    

###################
### plot 3, PIP ###
###################
    

dap2 <- dap%>%filter(gene==ens)
mapDF <- str_split(dap2$id_b37, ":", simplify=T)
dap2$chr <- as.numeric(mapDF[,1])
dap2$pos <- as.numeric(mapDF[,2])
dap2$pos2 <- (dap2$pos)/1e+03    
dap2$chr_pos_grch37 <- paste(mapDF[,1], mapDF[,2], sep="_")
dap2$rs <- gsub(".*;", "", mapDF[,4])
    
dap3 <- dap2%>%filter(chr==chr_i, pos>(pos_i-dtss1), pos<(pos_i+dtss2))%>%
    left_join(ldDF2, by="chr_pos_grch37")
    
dap3 <- dap3%>%mutate(gr2=ifelse(chr_pos_grch37==chr_pos, 0, gr))


##
anno <- dap3%>%filter(chr_pos_grch37==chr_pos)

 
p3 <- ggplot(dap3, aes(x=pos2, y=PIP))+
   geom_point(aes(colour=factor(gr2), size=factor(gr2),shape=factor(gr2)))+
   scale_size_manual("",
       values=c("0"=3.5, "1"=2,"2"=1.5, "3"=1.5, "4"=1.5, "5"=1.5, "6"=1.5), guide="none")+
   scale_color_manual("",
       values=c("1"="#D43F3A", "2"="#EEA236", "3"="#5CB85C", "4"="#46B8DA",
                "5"="#357EBD", "6"="#B8B8B8", "0"="#9632B8"),
       breaks=c(1, 2, 3, 4, 5),
       labels=c("1"="0.8~1", "2"="0.6~0.8", "3"="0.4~0.6",
                "4"="0.2~0.4", "5"="0~0.2"),
       guide=guide_legend(override.aes=list(size=2)))+
   scale_shape_manual("",
       values=c("0"=18, "1"=19, "2"=19, "3"=19, "4"=19, "5"=19, "6"=19), guide="none")+   
   geom_text_repel(data=anno, aes(x=pos2, y=PIP, label=rs),
                   color="#9632B8",
                   min.segment.length=0, segment.curvature=1e-20,
                   size=3.5, box.padding=0.8)+
   xlab(bquote("position (kb)"~"(chr"~.(chr_i)~")"))+
   ylab(bquote("Fine-map"~"("~italic(PIP)~")"))+ylim(0,1)+
   ## ggtitle(bquote(~italic(.(symbol))))+ 
   theme(plot.title=element_text(hjust=0.5, size=14),
      axis.text=element_text(size=10),
      axis.title=element_text(size=10),
      plot.margin=unit(c(5.5, 12, 5.5, 5.5), "pt"),
      legend.position=c(0.85,0.8),
      legend.key=element_blank(),
      legend.key.size=grid::unit(0.4, "cm"),
      legend.background=element_blank(),
      legend.box.background=element_blank(),
      legend.text=element_text(size=8),
      panel.background=element_blank(),
       panel.border=element_rect(fill=NA,color="black"),
       panel.grid.major=element_blank(),
       panel.grid.minor=element_blank())


    
##################
### gene track ###
##################
    
pos_min <- as.integer(pos_i-dtss1)
pos_max <- as.integer(pos_i+dtss2)
region2 <- paste(chr_i, pos_min, pos_max, sep="-")    
    
p4 <- AnnotationPlot(atac, region=region2)&
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())

ii <- length(p4$layers)
p4$layers[[ii]]$aes_params$fontface <- "italic"
p4$layers[[ii]]$aes_params$size <- 4
    
    
###
comb <- plot_grid(p1, p2, p3, p4, align="hv", axis="lr", nrow=4, rel_heights=c(0.6, 0.6, 0.6, 0.5))
figfn <- paste(outdir2, "Figure_chr17_", i, "_", ens, "_", symbol, "_1.comb.png", sep="")       
png(figfn, width=800, height=1050, res=120)
print(comb)
dev.off()

}
 
###END plot-1
 

###
###




###################
### plot-2
###################


contrast_ls <- list("LPS"=c("LPS", "CTRL"), "LPS-DEX"=c("LPS-DEX", "LPS"),
     "PHA"=c("PHA", "CTRL"), "PHA-DEX"=c("PHA-DEX", "PHA"))               
col_treat <- c("LPS"="#fb9a99", "LPS-DEX"="#e31a1c", "PHA"="#a6cee3", "PHA-DEX"="#1f78b4")
MCls_name <- c("Bcell"="B cell", "Monocyte"="Monocyte", "NKcell"="NK cell", "Tcell"="T cell")



trait <- sort(read.table("traits.txt")$V1)[3]
fn <- paste("./4_INTACT.outs/", trait, "_aloft_topPIP_union_combinfor.txt", sep="")
res <- read.table(fn, header=T, sep="\t")

x <- str_split(res$genetic_variant, ":", simplify=T)
res <- res%>%mutate(chr_pos_grch37=paste(x[,1], x[,2], sep="_"))

res3 <- res%>%dplyr::filter(FDR<0.1, FDR_twas<0.1, peaking_d==3, PIP>0.1)%>%arrange(FDR_twas)

isel <- c(1, 50)
conditions <- c("Bcell_LPS-DEX", "Monocyte_PHA")
tmp <- data.frame(irow=isel, "Gene"=res3$Gene[isel],
    "symbol"=res3$symbol[isel], "conditions"=conditions,
    "peaks"=c("5-141487761-141490047", "1-2486002-2488823"), ylims=c(100, 150),
    dtss_left=c(2500, 1000), dtss_right=c(1000,1000))

## ##
## opfn <- "./5_pub.outs/3_example_plots/geneList_param.txt"
## write.table(tmp, file=opfn, quote=F, row.names=F, sep="\t")

## tmp <- read.table("./5_pub.outs/3_example_plots/geneList_param.txt", header=T, sep="\t")

## dtss1 <- 5e+04
## dtss2 <- 5e+04
## i <- c(7, 11, 12, 14, 20, 21, 25)  
## for (i in  1:nrow(tmp)){
###
k <- 2
i <- tmp$irow[k]
dtss1 <- tmp$dtss_left[k]
dtss2 <- tmp$dtss_right[k]

ens <- res3$Gene[i]
symbol2 <- res3$symbol[i]
chr_pos <- res3$chr_pos_grch37[i]
chr_i <- as.numeric(gsub("_.*", "", chr_pos))
pos_i <- as.numeric(gsub(".*_", "",  chr_pos))
snp_i <- res3$genetic_variant[i]
    
cat(i, ens, symbol2, "\n")

###if (grepl("HLA", symbol)) next

#############################
### 0, calculate LD files ###    
#############################

###
### vcf file
genfn <- paste("./5_pub.outs/3_example_plots/vcf/", ens, "_", symbol2, "_", chr_i, "_dosages.txt", sep="")
gen <- read.table(genfn, header=F)
 
## ###
## gen2 <- gen%>%filter(V2>pos_i-dtss1, V2<pos_i+dtss2)
## isnp <- which(gen2$V2==pos_i)
 
### LD data frame    
## ldDF <- NULL
## g0 <- as.numeric(gen2[isnp,4:ncol(gen2)])    
## for (k in 1:nrow(gen2)){
##    ##
##    g1 <- as.numeric(gen2[k,4:ncol(gen2)])
##    r2 <- cor(g0, g1)^2
##    ## 
##    DF <- data.frame(chr=gen2[k,1], pos=gen2[k,2], genetic_variant=gen2[k,3], r2=round(r2, digits=3))
##    ldDF <- rbind(ldDF, DF) 
## }    

## ldDF2 <- ldDF%>%
##     mutate(chr_pos_grch37=paste(chr, pos, sep="_"))%>%
##     dplyr::select(chr_pos_grch37, r2)%>%
##     mutate(gr=case_when(r2>=0.8~"1", r2>=0.6&r2<0.8~"2", r2>=0.4&r2<0.6~"3",
##                         r2>=0.2&r2<0.4~"4", TRUE~"5"))



#################
### PIP plots ###
#################
    
 
## dap2 <- dap%>%filter(gene==ens)
## mapDF <- str_split(dap2$id_b37, ":", simplify=T)
## dap2$chr <- as.numeric(mapDF[,1])
## dap2$pos <- as.numeric(mapDF[,2])
## dap2$chr_pos_grch37 <- paste(mapDF[,1], mapDF[,2], sep="_")
## dap2$rs <- gsub(".*;", "", mapDF[,4])

## dap3 <- dap2%>%filter(chr==chr_i, pos>(pos_i-dtss1), pos<(pos_i+dtss2))%>%
##     mutate(pos2=pos/1e+03)%>%
##     left_join(ldDF2, by="chr_pos_grch37")
    
## dap3 <- dap3%>%mutate(gr2=ifelse(is.na(r2), 6, gr))%>%mutate(gr2=ifelse(r2==1, 0, gr2))  

## ##
## anno <- dap3%>%filter(chr_pos_grch37==chr_pos)

 
## pp1 <- ggplot(dap3, aes(x=pos2, y=PIP))+
##    geom_point(aes(colour=factor(gr2), size=factor(gr2),shape=factor(gr2)))+
##    scale_size_manual("",
##        values=c("0"=3, "1"=1.5,"2"=2, "3"=2, "4"=2, "5"=2, "6"=2), guide="none")+
##    scale_color_manual("",
##        values=c("1"="#D43F3A", "2"="#EEA236", "3"="#5CB85C", "4"="#46B8DA",
##                 "5"="#357EBD", "6"="#B8B8B8", "0"="#9632B8"),
##        breaks=c(1, 2, 3, 4, 5),
##        labels=c("1"="0.8~1", "2"="0.6~0.8", "3"="0.4~0.6",
##                 "4"="0.2~0.4", "5"="0~0.2"),
##        guide=guide_legend(override.aes=list(size=2)))+
##    scale_shape_manual("",
##        values=c("0"=18, "1"=19, "2"=19, "3"=19, "4"=19, "5"=19, "6"=19), guide="none")+   
##    geom_text_repel(data=anno, aes(x=pos2, y=PIP, label=rs),
##                    color="#9632B8", size=3.5, box.padding=0.8)+
##    xlab(bquote("position (kb)"~"(chr"~.(chr_i)~")"))+
##    ylab(bquote("Fine-map"~"("~italic(PIP)~")"))+ylim(0,1)+
##    ggtitle(bquote(~italic(.(symbol))))+ 
##    theme(plot.title=element_blank(), ##(hjust=0.5, size=14),
##       axis.text=element_text(size=10),
##       axis.title=element_text(size=10),      
##       plot.margin=unit(c(5.5, 12, 5.5, 5.5), "pt"),
##       legend.position=c(0.85,0.7),
##       legend.key=element_blank(),
##       legend.key.size=grid::unit(0.4, "cm"),
##       legend.background=element_blank(),
##       legend.box.background=element_blank(),
##       legend.text=element_text(size=8),     
##       panel.background=element_blank(),
##        panel.border=element_rect(fill=NA,color="black"),
##        panel.grid.major=element_blank(),
##        panel.grid.minor=element_blank()) 



#############################################
### 2, Show specific genomic track region ###
#############################################
 

### ATAC, Chromatin accessibility
## default 1000
pos_min <- pos_i-dtss1
pos_min <- ifelse(pos_min<0, 0, as.integer(pos_min))
pos_max <- as.integer(pos_i+dtss2)
 

##condition <- DFanno%>%filter(symbol==symbol2)%>%pull(conditions)%>%unique()
oneMCl <- gsub("_.*", "", tmp$conditions[k])
oneMCl2 <- MCls_name[oneMCl]

ii <- gsub(".*_", "", tmp$condition[k])
i1 <- contrast_ls[[ii]][1]
i2 <- contrast_ls[[ii]][2]
    
 
atac2 <- subset(atac, MCls==oneMCl&treat%in%contrast_ls[[ii]])
 
peak_i <- tmp$peaks[k]
peak_i <- StringToGRanges(peak_i)
    
region <- paste(chr_i, pos_min, pos_max, sep="-") 
pp2 <- CoveragePlot(atac2, region=region, show.bulk=F, peaks=F, annotation=F, tile=F, 
   group.by="treat", region.highlight=peak_i, extend.upstream=0e+00, extend.downstream=0e+00, links=F,
   heights=c(1, 0.8))&
   ## scale_fill_manual(values=c("Bcell"="#4daf4a", "Monocyte"="#984ea3",
   ##     "NKcell"="#aa4b56", "Tcell"="#ffaa00","DC"="#828282"))&
   scale_fill_manual(values=c("LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
                              "PHA"="#a6cee3", "PHA-DEX"="#1f78b4"))& 
   annotate("point",x=pos_i, y=1, shape=17, size=3.5, color="blue")&
   xlab(bquote(~"position (bp)"~"(chr"~.(chr_i)~")"))&
   ggtitle(oneMCl2)&
   theme(legend.position="none",
         plot.title=element_text(hjust=0.5, size=12),
         axis.title=element_text(size=10),
         ##axis.text.x=element_blank(),
         ##axis.ticks.x=element_blank(),
         axis.ticks.y=element_blank(),
         strip.text.y.left=element_blank())


##############################################################
### genomic track region for different genotype separately ###
##############################################################
    
id6 <- unique(atac$SNG.BEST.GUESS)
ylim2 <- tmp$ylims[k]

###
### genotype
str2 <- str_split(gsub(";.*", "", snp_i), ":", simplify=T)
ref <- str2[1,3]
alt <- str2[1,4]
genotype <- c(paste0(ref, ref), paste0(ref, alt), paste0(alt, alt)) 
names(genotype) <- c("0", "1", "2")
###
 
gen_i <- gen%>%filter(V3==snp_i)
genDF <- data.frame(ID=sampleID, gr=as.character(round(gen_i[1,4:254])))%>%
    filter(ID%in%id6)%>%
    mutate(genotype2=genotype[gr])


### ATAC, Chromatin accessibility
##condition <- DFanno%>%filter(symbol==symbol2)%>%pull(conditions)%>%unique()
## oneMCl <- "Monocyte" #gsub("_.*", "", condition)
## ii <- "LPS-DEX"  ## gsub(".*_", "", condition)
## i1 <- contrast_ls[[ii]][1]
## i2 <- contrast_ls[[ii]][2]

 
### atac data
atac2 <- subset(atac, MCls==oneMCl&treat%in%contrast_ls[[ii]])
meta <- atac2@meta.data
meta2 <- meta%>%left_join(genDF, by=c("SNG.BEST.GUESS"="ID"))
rownames(meta2) <- meta2$NEW_BARCODE
atac2 <- AddMetaData(atac2, meta2)

 
### ref-00, hete-1 and alt-2
pp_ls <- lapply(0:2, function(kk){
  ##
  gg <- genotype[as.character(kk)]
  atac3 <- try(subset(atac2, genotype2==gg), silent=T)

  ##  
  if ( class(atac3)!="try-error"){
###      
   nindi <- length(unique(atac3$SNG.BEST.GUESS))     
###
  pos_min <- pos_i-dtss1
  pos_min <- ifelse(pos_min<0, 0, as.integer(pos_min))
  pos_max <- as.integer(pos_i+dtss2)
    
  region <- paste(chr_i, pos_min, pos_max, sep="-")
    
    
  p0 <- CoveragePlot(atac3, region=region, show.bulk=F, peaks=F, annotation=F, tile=F, 
     group.by="treat", ymax=ylim2, extend.upstream=0e+00, extend.downstream=0e+00, links=F,
     heights=c(1, 0.8), region.highlight=peak_i)&
     ## scale_fill_manual(values=c("Bcell"="#4daf4a", "Monocyte"="#984ea3",
     ##     "NKcell"="#aa4b56", "Tcell"="#ffaa00","DC"="#828282"))&
     scale_fill_manual(values=c("LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
                                "PHA"="#a6cee3", "PHA-DEX"="#1f78b4"))& 
     annotate("point",x=pos_i, y=1, shape=17, size=3.5, color="blue")&
     xlab(bquote(~"position (bp)"~"(chr"~.(chr_i)~")"))&
     ggtitle(bquote(.(oneMCl2)~"("~.(gg)~","~.(nindi)~")"))&
     theme(legend.position="none",
           plot.title=element_text(hjust=0.5, size=12),
           axis.title.x=element_text(size=10),
           axis.title.y=element_text(size=8),
           ##axis.text.x=element_blank(),
           ##axis.ticks.x=element_blank(),
           axis.ticks.y=element_blank(),
           strip.text.y.left=element_blank())
    ###
    }else{
       p0 <- NULL
    }   
    p0      
})


sel <- !sapply(pp_ls, is.null)

    

###     
### output plots

if ( sum(sel)==3){
### 3-genotype    
comb <- plot_grid(pp2, pp_ls[[1]], pp_ls[[2]], pp_ls[[3]],
    align="hv", axis="lr", nrow=4, rel_heights=c(1, 1, 1, 1))
  
figfn <- paste(outdir2, "Figure_", i, "_", ens, "_", symbol2, "_2.peak.png", sep="")       
png(figfn, width=800, height=800, res=120)
print(comb)
dev.off()
###
}else{    
### 2-genotype
pi_1 <- which(sel)[1]
pi_2 <- which(sel)[2]    
comb <- plot_grid(pp2, pp_ls[[pi_1]], pp_ls[[pi_2]], 
    align="hv", axis="lr", nrow=3, rel_heights=c(1, 1, 1))
  
figfn <- paste(outdir2, "Figure_", i, "_", ens, "_", symbol2, "_2.peak.png", sep="")       
png(figfn, width=800, height=720, res=120)
print(comb)
dev.off()
    
} ##

    

    
  
## pos_min <- as.integer(pos_i-dtss1)
## pos_max <- as.integer(pos_i+dtss2)
## region2 <- paste(chr_i, pos_min, pos_max, sep="-")    
    
## p5 <- AnnotationPlot(atac2, region=region2)&
##     theme(axis.title.x=element_blank(),
##           axis.text.x=element_blank(),
##           axis.ticks.x=element_blank())
    
## ## figfn <- paste(outdir2, "Figure", i, ".4_", ens, "_", symbol, "_peaks.png", sep="")       
## ## png(figfn, width=600, height=400, res=120)
## ## print(p4)
## ## dev.off()

## ###
## ###
## comb <- plot_grid(p1, p2, p3, p4, p5, align="hv", axis="lr", nrow=5, rel_heights=c(0.7, 0.6, 0.6, 0.6, 0.4))
## figfn <- paste(outdir2, "Figure", i, ".1_", ens, "_", symbol, "_comb.png", sep="")       
## png(figfn, width=800, height=1250, res=120)
## print(comb)
## dev.off()

## }


######################
### SNP annotation ###
######################


### SNP hit by motifs

prefix <- "/nfs/rprdata/julong/sc-atac/genetic.analysis_torus_2021-10-14/SNPannotation/"
fn <- paste(prefix, "motifList2020.txt", sep="")
motifList <- read.table(fn)$V1
snpmotif <- lapply(motifList, function(motif){
###
## cat(motif, "\n") 
   fn <- paste(prefix, "annot_jaspar2020/allsnp_", motif, ".bed.gz", sep="")   
   xx <- try(fread(fn, header=F, data.table=F, stringsAsFactors=F), silent=T)
   if ( class(xx)!="try-error"){ 
      xx$motifs <- motif
   }else{
      xx <- NA
   }   
   xx
})
snpmotif <- do.call(rbind, snpmotif[!is.na(snpmotif)])
snpmotif <- snpmotif%>%mutate(chr_pos=paste(V1, V2, sep="_"))


###
### response motifs
## contrast <- c("LPS", "LPS-DEX", "PHA", "PHA-DEX")
## MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
## ###
## combs <- paste(rep(MCls, each=4), rep(contrast, times=4), sep="_")
## ##
## motifRsps <- map_dfr(combs, function(ii){
##    ##
##    prefix <- "/nfs/rprdata/julong/sc-atac/analyses.2021-02-05/3_motif/2_motif.activities.outs/MotifList/" 
##    fn <- paste(prefix, "1_comb_", ii, "_motif.txt", sep="")
##    tmp <- read.table(fn)
##    names(tmp) <- "motif.id"
##    tmp$conditions <- ii
##    tmp
## })    


###
### differential motif
fn <- "/nfs/rprdata/julong/sc-atac/analyses.2021-02-05/3_motif/2_motif.activities.outs/3_motif.diff.results.rds"
resDiff <- read_rds(fn)%>%
    mutate(conditions=paste(MCls, contrast, sep="_"), comb=paste(gene, MCls, contrast, sep="_"))
resDiff2 <- resDiff%>%dplyr::select(conditions, gene, motif, beta, pval, qval)



res3 <- res2%>%filter(peaking_d==3, PIP>0.7|FDR_intact<0.1)

####
#### SNP annotation table
DFanno <- NULL
for (i in 1:nrow(res3)){
  ##   
  chr_pos1 <- res3$chr_pos_grch37[i]
  symbol2 <- res3$symbol[i]
  snp <- res3$genetic_variant[i]
##  if ( grepl("HLA", symbol)) next
  cat(symbol2, "\n")  
  ##  
  motif0 <- snpmotif%>%filter(chr_pos==chr_pos1)%>%pull(motifs)
  ###DF <- motifMCls%>%filter(motif.id%in%motif0)
  ##
  ## DF2 <- motifRsps%>%filter(motif.id%in%motif0)%>%mutate(comb=paste(motif.id, conditions, sep="_"))    
  ## DF2 <- DF2%>%left_join(resDiff2, by="comb") ##%>%slice_max(order_by=abs(beta), n=1)
  DF2 <- resDiff2%>%filter(gene%in%motif0)%>%arrange(gene)  
     
  ###DF2 <- DF2%>%full_join(DF, by="motif.id")
  ##
  DF2 <- DF2%>%mutate(symbol=symbol2, genetic_variant=snp)
  DFanno <- rbind(DFanno, DF2)  
}    
  
## DFanno <- DFanno%>%drop_na(conditions)

opfn <- "./5_pub.outs/3_example_plots/ALOFT_union/zzz_motif_annot.xlsx"
write.xlsx(DFanno, file=opfn, overwrite=T)
  




##########################
### Differential peaks ###
##########################

fn <- "../analyses.2021-02-05/2_Differential/1.3_DiffPeak.outs/3.0_DESeq_indi.results.rds"
resDAR <- read_rds(fn)%>%mutate(is_sig=ifelse(p.adjusted<0.1&abs(estimate)>0.5, sign(estimate), 0))

##
gene <- unique(resDAR$gene)
cvt <- str_split(gene, "-", simplify=T) 
DF <- data.frame(gene=gene, chr=as.numeric(cvt[,1]), pos1=as.integer(cvt[,2]), pos2=as.integer(cvt[,3]))

##
resAnno <- NULL
for (i in 1:nrow(res3)){
   ##
   chr_pos <- res3$chr_pos_grch37[i]
   chri <- as.numeric(gsub("_.*", "", chr_pos))
   posi <- as.integer(gsub(".*_", "", chr_pos))
   symbol_i <- res3$symbol[i]
   SNP_i <- res3$genetic_variant[i]
   PIP_i <- res3$PIP[i]
   pval_gwas_i <- res3$pval_gwas[i]
    
   ###
   peakSel <- DF%>%filter(chr==chri, pos1<=posi, pos2>=posi)%>%pull(gene)
   ##
   resDAR2 <- resDAR%>%filter(gene%in%peakSel)%>%dplyr::select(-baseMean, -statistic)
   resDAR2 <- resDAR2%>%
       mutate(symbol=symbol_i, genetic_variant=SNP_i, PIP=PIP_i, chr_pos_grch37=chr_pos, pval_gwas=pval_gwas_i)
   ##
   resAnno <- rbind(resAnno, resDAR2)
}

opfn <- "./5_pub.outs/3_example_plots/ALOFT_union/zzz_DAR_annot.xlsx"
write.xlsx(resAnno, file=opfn, overwrite=T)


###
## resAnno2 <- resAnno%>%mutate(conditions=paste(MCls, contrast, sep="_"))%>%
##     group_by(symbol)%>%slice_min(order_by=p.value, n=4)%>%ungroup()%>%arrange(pval_gwas)
 
## opfn2 <- "./4.3_example.outs/ALOFT_union/genes_list3/zzz2_DAR_short.xlsx"
## openxlsx::write.xlsx(resAnno2, file=opfn2, overwrite=T)





#########################################
### for speific motifs used for plots ###
#########################################

for (symbol0 in geneSel){    
##symbol0 <- "STIM2"
resMotif <- read.xlsx("./4.3_example.outs/ALOFT_union/zzz_motif_annot.xlsx")

x <- resMotif%>%filter(symbol==symbol0)
x <- x%>%arrange(desc(abs(beta)))%>%
    mutate(beta=round(beta,3))%>% ## pval=round(pval,3), qval=round(qval, 3))%>%
    dplyr::select(conditions, motif, beta, pval, qval)
###
opfn <- paste("5_pub.outs/3_example_plots/ALOFT_union/tmp_", symbol0, "_motif.xlsx", sep="")
write.xlsx(x, opfn, overwrite=T)


###
### 
fn <- "./4.3_example.outs/ALOFT_union/zzz_DAR_annot.xlsx"
resDAR <- openxlsx::read.xlsx(fn)
tmp2 <- resDAR%>%filter(symbol==symbol0)%>%arrange(desc(abs(estimate)))%>%
    mutate(conditions=paste(MCls, contrast, sep="_"),
           LFC=round(estimate, 3), stderror=round(stderror, 3),
           pval=p.value, qval=p.adjusted)%>%
           ##pval=round(p.value, 3), qval=round(p.adjusted, 3))%>%
    filter(MCls!="DC")%>%
    dplyr::select(conditions, gene, LFC, stderror, pval, qval, is_sig)    

####
opfn <- paste("5_pub.outs/3_example_plots/ALOFT_union/tmp_", symbol0, "_DAR.xlsx", sep="")
write.xlsx(tmp2, opfn, overwrite=T)
}


### motif
## fn <- "./4.3_example.outs/ALOFT_union/genes_list3/zzz_annot.xlsx"
## resmotif <- read.xlsx(fn)
## x2 <- resmotif%>%filter(motif=="PAX5", symbol==symbol_i)%>%arrange(pval)%>%
##     mutate(beta=round(beta,3), pval=round(pval,3), qval=round(qval, 3))%>%
##     dplyr::select(conditions, motif, beta, pval, qval, symbol)
## ##
## opfn <- paste(outdir2, "tmp_", symbol_i, "_motif.xlsx", sep="")
## write.xlsx(x2, opfn, overwrite=T)



## fn <- "/nfs/rprdata/julong/sc-atac/genetic.analysis_torus_2021-10-14/SNPannotation/liftOver/SCAIP_final_bed.gz"
## x <- fread(fn, header=F)


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
