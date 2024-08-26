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
 

###                            
### candidate genes 
fn <- "./5_pub.outs/2_supp_tables/TableS5_1_asthma-risk-genes_ALOFT.txt.gz"
res <- read.table(fn, header=T, sep="\t")
 

res2 <- res%>%filter(PIP>0.1, FDR_intact<0.1, FDR_twas<0.1)

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
### gwas

trait <- sort(read.table("traits.txt")$V1)[3]
fn <- paste("./gwas_impute/", trait, "_impute.txt.gz", sep="")
gwas <- fread(fn, header=T, data.table=F)
gwas <- gwas%>%mutate(chr_pos_grch37=pos37[id_b38_0])

gwas <- gwas%>%dplyr::select(chr_pos_grch38=id_b38_0, chr, zscore, pval, chr_pos_grch37)
 

###
### dap results 
fn <- "/nfs/rprdata/julong/sc-atac/twas_analysis_2022-10-16/ALOFT_results/aloft_allSNPs_union.txt.gz"

dap <- fread(fn, header=F, data.table=F)
names(dap) <- c("gene", "genetic_variant", "PIP")

###
### FastQTL results
fn <- "/nfs/rprdata/julong/sc-atac/twas_analysis_2022-10-16/ALOFT_results/PC1-18.nominals.eQTL.txt.gz"
fast <- fread(fn, data.table=F)
names(fast) <- c("gene", "genetic_variant", "DTSS", "pval", "beta")



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


outdir2 <- "./5_pub.outs/3_example_plots/Example2_aloft/"
if (!file.exists(outdir2)) dir.create(outdir2, showWarnings=F, recursive=T)


###
###
fn <- "./5_pub.outs/2_supp_tables/TableS5_1_asthma-risk-genes_ALOFT.txt.gz"
res <- fread(fn, header=T, sep="\t", data.table=F)

## main
res2 <- res%>%dplyr::filter(FDR_intact<0.1, FDR_twas<0.1, PIP>0.1, peaking_d==3, pval_twas<1e-07)
 
## opfn <- paste(outdir, "geneList.xlsx", sep="")
## write.xlsx(res2, file=opfn)
## opfn <- paste(outdir, "geneList.txt", sep="")
## write.table(res2, file=opfn, quote=F, row.names=F)

## geneSel <- c("KIF1B", "TNFRSF14", "SIK2", "NDFIP1",  "TRAT1", "HAUS8", "MPHOSPH9")[c(1,2,4)]
## res3 <- res3%>%filter(symbol%in%geneSel)
## res3 <- res3[c(3,4,8),]
    ##filter(peaking_d==dd, PIP>0.5, FDR_intact<0.1)


res3 <- res2
x <- str_split(res3$genetic_variant, ":", simplify=T)
res3$chr_pos_grch37 <- paste(x[,1], x[,2], sep="_")


###############
### plots-1 ###
###############
 
##dtss <- 2.5e+04 ## default 5e+04, 10e+04
dtss1 <- 5e+04
dtss2 <- 5e+04

 
for (i in  1:nrow(res3)){
###    

 
ens <- res3$Gene[i]
symbol2 <- res3$symbol[i]

chr_pos <- res3$chr_pos_grch37[i]
chr_i <- as.numeric(gsub("_.*", "", chr_pos))
pos_i <- as.numeric(gsub(".*_", "",  chr_pos))
snp_i <- res3$genetic_variant[i]
dd <- res3$peaking_d[i]
    
cat(i, ens, symbol2, "\n")

###if (grepl("HLA", symbol)) next

#########################
### local zoom  plots ###
#########################

    
#############################
### 0, calculate LD files ###    
#############################

###
### vcf file
genfn <- paste("./5_pub.outs/3_example_plots/vcf/", ens, "_", symbol2, "_", chr_i, "_dosages.txt", sep="")
gen <- read.table(genfn, header=F)

###
gen2 <- gen%>%filter(V2>(pos_i-dtss1), V2<(pos_i+dtss2))
isnp <- which(gen2$V3==snp_i)

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
    dplyr::select(genetic_variant, chr_pos_grch37, r2)%>%
    mutate(gr=case_when(r2>=0.8~"1", r2>=0.6&r2<0.8~"2", r2>=0.4&r2<0.6~"3",
                        r2>=0.2&r2<0.4~"4", TRUE~"5"))

ldDF2 <- ldDF2%>%mutate(rs=gsub(".*;", "", genetic_variant))

    

    
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
fast2$chr <- as.integer(mapDF[,1])
fast2$pos <- as.integer(mapDF[,2])
       
fast3 <- fast2%>%filter(chr==chr_i, pos>(pos_i-dtss1), pos<(pos_i+dtss2))%>%
    left_join(ldDF2, by="genetic_variant")
    
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
mapDF <- str_split(dap2$genetic_variant, ":", simplify=T)
dap2$chr <- as.integer(mapDF[,1])
dap2$pos <- as.integer(mapDF[,2])
dap2$pos2 <- dap2$pos/1e+03
    
    
dap3 <- dap2%>%filter(chr==chr_i, pos>(pos_i-dtss1), pos<(pos_i+dtss2))%>%
    left_join(ldDF2, by="genetic_variant")
    
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
if (ii>0){    
p4$layers[[ii]]$aes_params$fontface <- "italic"
p4$layers[[ii]]$aes_params$size <- 4
}    
    
###
comb <- plot_grid(p1, p2, p3, p4, align="hv", axis="lr", nrow=4, rel_heights=c(0.6, 0.6, 0.6, 0.5))
figfn <- paste(outdir2, "Figure_anno", dd, "_", i, "_", ens, "_", symbol2, "_", chr_i, "_1.comb.png", sep="")       
png(figfn, width=800, height=1050, res=120)
print(comb)
dev.off()

}


###END plot-1
###





###################
### plot-2
###################


outdir2 <- "./5_pub.outs/3_example_plots/Example2_aloft_atac/"
if (!file.exists(outdir2)) dir.create(outdir2, showWarnings=F, recursive=T)


contrast_ls <- list("LPS"=c("LPS", "CTRL"), "LPS-DEX"=c("LPS-DEX", "LPS"),
     "PHA"=c("PHA", "CTRL"), "PHA-DEX"=c("PHA-DEX", "PHA"))               
col_treat <- c("LPS"="#fb9a99", "LPS-DEX"="#e31a1c", "PHA"="#a6cee3", "PHA-DEX"="#1f78b4")
MCls_name <- c("Bcell"="B cell", "Monocyte"="Monocyte", "NKcell"="NK cell", "Tcell"="T cell")


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






fn <- "./5_pub.outs/2_supp_tables/TableS5_2_asthma-risk-genes_in_response_ALOFT.txt.gz"
res <- fread(fn, header=T, sep="\t", data.table=F)
res <- res%>%filter(PIP>0.1)

## x <- str_split(res$genetic_variant, ":", simplify=T)
## res <- res%>%mutate(chr_pos_grch37=paste(x[,1], x[,2], sep="_"))
## res3 <- res%>%dplyr::filter(FDR<0.1, FDR_twas<0.1, peaking_d==3, PIP>0.1)%>%arrange(FDR_twas)


isel <- c(1, 4, 7, 8, 10, 45)

conditions <- c("Bcell_LPS-DEX", "Tcell_PHA-DEX", "Tcell_LPS-DEX",
                "Monocyte_PHA", "Tcell_PHA", "Monocyte_PHA")



plotDF <- data.frame(irow=isel, "Gene"=res$Gene[isel], "symbol"=res$symbol[isel], "conditions"=conditions,
   "peaks"=res$peakAnno[isel],"chr_pos_grch37"=res$chr_pos_grch37[isel],"genetic_variant"=res$genetic_variant[isel])
    
## , ylims=c(100, 150),
##     dtss_left=c(2500, 1000), dtss_right=c(1000,1000))

## ##
## opfn <- "./5_pub.outs/3_example_plots/geneList_param.txt"
## write.table(tmp, file=opfn, quote=F, row.names=F, sep="\t")

## tmp <- read.table("./5_pub.outs/3_example_plots/geneList_param.txt", header=T, sep="\t")

## dtss1 <- 5e+04
## dtss2 <- 5e+04
## i <- c(7, 11, 12, 14, 20, 21, 25)  
## for (i in  1:nrow(tmp)){
###


for (k in 1:nrow(plotDF)){

###    

i <- plotDF$irow[k]


ens <- plotDF$Gene[k]
symbol2 <- plotDF$symbol[k]
chr_pos <- plotDF$chr_pos_grch37[k]
chr_i <- as.numeric(gsub("_.*", "", chr_pos))
pos_i <- as.numeric(gsub(".*_", "",  chr_pos))
snp_i <- plotDF$genetic_variant[k]


    
peak0 <- plotDF$peaks[k]

peak_sp <- unlist(str_split(peak0, "-"))
pos_min <- as.integer(peak_sp[2])-500
pos_max <- as.integer(peak_sp[3])+500 
    
## dtss1 <- as.integer(pos_i)-as.integer(peak_sp[2])+1000  ##tmp$dtss_left[k]
## dtss2 <- as.integer(pos_i)-as.integer(peak_sp[3])+100
 

    
cat(i, ens, symbol2, "\n")

###if (grepl("HLA", symbol)) next

#############################
### 0, calculate LD files ###    
#############################

###
### vcf file
genfn <- paste("./5_pub.outs/3_example_plots/vcf/", ens, "_", symbol2, "_", chr_i, "_dosages.txt", sep="")
gen <- read.table(genfn, header=F)

    
id6 <- unique(atac$SNG.BEST.GUESS)

###
### genotype
str2 <- unlist(str_split(gsub(";.*", "", snp_i), ":"))
ref <- str2[3]
alt <- str2[4]
genotype <- c(paste0(ref, ref), paste0(ref, alt), paste0(alt, alt)) 
names(genotype) <- c("0", "1", "2")
###
 
gen_i <- gen%>%filter(V3==snp_i)
genDF <- data.frame(ID=sampleID, gr=as.character(round(gen_i[1,4:254])))%>%
    filter(ID%in%id6)%>%
    mutate(genotype2=genotype[gr])


    

    
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

## pos_min <- pos_i-dtss1
## pos_min <- ifelse(pos_min<0, 0, as.integer(pos_min))
## pos_max <- as.integer(pos_i+dtss2)
 

##condition <- DFanno%>%filter(symbol==symbol2)%>%pull(conditions)%>%unique()
oneMCl <- gsub("_.*", "", plotDF$conditions[k])
oneMCl2 <- MCls_name[oneMCl]

ii <- gsub(".*_", "", plotDF$condition[k])
i1 <- contrast_ls[[ii]][1]
i2 <- contrast_ls[[ii]][2]
     
atac2 <- subset(atac, MCls==oneMCl&treat%in%contrast_ls[[ii]])

meta <- atac2@meta.data
meta2 <- meta%>%left_join(genDF, by=c("SNG.BEST.GUESS"="ID"))
rownames(meta2) <- meta2$NEW_BARCODE
atac2 <- AddMetaData(atac2, meta2)    

    
peak_i <- plotDF$peaks[k]
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

    
### ATAC, Chromatin accessibility
##condition <- DFanno%>%filter(symbol==symbol2)%>%pull(conditions)%>%unique()
## oneMCl <- "Monocyte" #gsub("_.*", "", condition)
## ii <- "LPS-DEX"  ## gsub(".*_", "", condition)
## i1 <- contrast_ls[[ii]][1]
## i2 <- contrast_ls[[ii]][2]



 
### ref-00, hete-1 and alt-2
pp_ls <- lapply(0:2, function(kk){
  ##
  gg <- genotype[as.character(kk)]
  atac3 <- try(subset(atac2, genotype2==gg), silent=T)

  ##  
  if ( class(atac3)!="try-error"){
###      
   nindi <- length(unique(atac3$SNG.BEST.GUESS))     

    

 ## ymax=ylim2,
  region <- paste(chr_i, pos_min, pos_max, sep="-")    
  p0 <- CoveragePlot(atac3, region=region, show.bulk=F, peaks=F, annotation=F, tile=F, 
     group.by="treat",  extend.upstream=0e+00, extend.downstream=0e+00, links=F,
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
  
figfn <- paste(outdir2, "Figure_anno3_", i, "_", ens, "_", symbol2, ".peak.png", sep="")       
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
  
figfn <- paste(outdir2, "Figure_anno3_", i, "_", ens, "_", symbol2, ".peak.png", sep="")       
png(figfn, width=800, height=720, res=120)
print(comb)
dev.off()
    
} ##

}  ### End gene-loop     

    

    
  
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
### Summary tables ###
######################

outdir2 <- "./5_pub.outs/3_example_plots/Example2_infor/"
if (!file.exists(outdir2)) dir.create(outdir2, showWarnings=F, recursive=T)


fn <- "./5_pub.outs/2_supp_tables/TableS5_1_asthma-risk-genes_ALOFT.txt.gz"
res <- fread(fn, header=T, sep="\t", data.table=F)

res2 <- res%>%filter(pval_twas<1e-07)

### select by Francesca and Roger
geneSel <- c("NDFIP1", "B3GALT6", "DEF6", "PTK7", "NCOR1", "TNFRSF14")

res2 <- res%>%dplyr::filter(FDR_intact<0.1, FDR_twas<0.1, PIP>0.1, peaking_d==3, symbol%in%geneSel)

### example genes
## res2 <- res%>%dplyr::filter(FDR_intact<0.1, FDR_twas<0.1, PIP>0.1, peaking_d==3, pval_twas<1e-07)
## res0 <- res%>%filter(symbol%in%c("MPHOSPH9", "TNFRSF14"))
## res2 <- rbind(res2, res0)

gene_val <- 1:nrow(res2)
names(gene_val) <- res2$Gene


###
### Differential accessbility

fn <- "./5_pub.outs/2_supp_tables/TableS5_3_asthma-risk-genes_peak.txt.gz"
resDP <- fread(fn, header=T, sep="\t", data.table=F)

resDP2 <- resDP%>%filter(Gene%in%res2$Gene, MCls!="DC")%>%
    mutate(condition=paste(MCls, contrast, sep="_"), gene_value=as.numeric(gene_val[Gene]),
           baseMean=round(baseMean, 3),
           estimate=round(estimate, 3), stderror=round(stderror, 3),
           pval=round(p.value, 3), qval=round(p.adjusted, 3))

resDP2 <- resDP2%>%
    dplyr::select(Gene, symbol, genetic_variant, condition, peak, baseMean,
       LFC=estimate,  stderror, pval, qval, gene_value)

resDP3 <- resDP2%>%       ###group_by(Gene)%>%slice_max(order_by=abs(LFC), n=10)%>%ungroup()%>%
    arrange(gene_value, qval)
 

opfn <- paste(outdir2, "Table_DP.infor.xlsx", sep="")
write.xlsx(resDP3, file=opfn, showNA=T)



###########################
### Differential motifs ###
###########################



### Binding score
fn <- "./5_pub.outs/2_supp_tables/TableS5_2_asthma-risk-genes_in_response_ALOFT.txt.gz"
x <- fread(fn, header=T, sep="\t")
x <- x%>%filter(Gene%in%res2$Gene)

###
cv <- NULL
for (i in 1:nrow(x)){
   ##
   mm0 <- unlist(str_split(x$motifAnnot[i], ";"))
   cv0 <- as.data.frame(str_split(mm0, "_", simplify=T))
   names(cv0) <- c("motif_ID", "score_ref", "score_alt") 
   cv0$Gene <- x$Gene[i]
   cv0$gene_motif <- paste(cv0$Gene, cv0$motif_ID, sep="_") 
   cv <- rbind(cv, cv0)
}    

th0 <- log2(20)

cv2 <- cv%>%mutate(diff=abs(as.numeric(score_ref)-as.numeric(score_alt)))%>%
       dplyr::select(gene_motif, diff) 
   ###group_by(Gene)%>%slice_max(order_by=diff, n=3)%>%ungroup()%>%


     


fn <- "./5_pub.outs/2_supp_tables/TableS5_4_asthma-risk-genes_motif.txt.gz"
resMotif <- fread(fn, header=T, sep="\t", data.table=F)

resMotif2 <- resMotif%>%filter(Gene%in%res2$Gene, MCls!="DC")%>%
    mutate(gene_value=as.numeric(gene_val[Gene]), motif_ID, motif_name, 
           beta=round(beta, 3), stderr=round(stderr, 3),
           pval=round(pval, 3), qval=round(qval, 3), gene_motif=paste(Gene, motif_ID, sep="_"))

## resMotif2 <- resMotif2%>%group_by(gene_motif)%>%slice_max(order_by=abs(beta), n=5)%>%ungroup()
 
resMotif2 <- resMotif2%>%inner_join(cv2, by="gene_motif")


resMotif2 <- resMotif2%>%
    dplyr::select(gene_motif, symbol, genetic_variant, condition=conditions, motif_name, 
       beta,  stderr, pval, qval, score_diff=diff, gene_value)

resMotif3 <- resMotif2%>%arrange(gene_value, desc(score_diff), qval)


####
fn2 <- paste(outdir2, "Table_DP.infor.xlsx", sep="")
DP <- read.xlsx(fn2)

resMotif_new <- NULL
for ( gene0 in unique(resMotif3$symbol)){
    ##
    conSel <- DP%>%filter(symbol==gene0, qval<0.1)%>%pull(condition)%>%unique()
    if ( length(conSel)==0){
       ###
        DP0 <- DP%>%filter(symbol==gene0)
        conSel <- DP0$condition[1]
    }
    
    res0 <- resMotif3%>%filter(symbol==gene0, condition%in%conSel)
    motifs <- unique(res0$motif_name)
    res0 <- map_dfr(motifs, function(ii){
        ##
        tmp <- res0%>%filter(motif_name==ii)
        rownames(tmp) <- tmp$condition
        tmp <- tmp[conSel,]
        tmp
    })    
    resMotif_new <- rbind(resMotif_new, res0)
}    

 
opfn <- paste(outdir2, "Table_Motif.infor.xlsx", sep="")
write.xlsx(resMotif_new, file=opfn, showNA=T)



###
### motif
fn <- "./5_pub.outs/3_example_plots/plots_pub_final/Table_motif.infor.txt"
x <- read.table(fn, header=T)

 

###################
### cluster PIP ###
###################

geneSel <- c("DEF6", "TNFRSF14", "NDFIP1")

fn <- "./5_pub.outs/2_supp_tables/TableS5_2_asthma-risk-genes_in_response_ALOFT.txt.gz"
res <- fread(fn, header=T, sep="\t")
res3 <- res%>%filter(symbol%in%geneSel)

## res3 <- res2%>%dplyr::filter(symbol%in%c("GSDMB", "ORMDL3", "IKZF3", "IL4", "FADS2"))

i <- 2   
ens <- res3$Gene[i]
symbol2 <- res3$symbol[i]

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


fn <- "./5_pub.outs/3_example_plots/plots_pub_final/Table_motif.infor.txt"
DF_motif <- read.table(fn, header=T, sep="\t")


mm <- unlist(str_split(res3$motifAnnot[i], ";"))
x <- str_split(mm, "_", simplify=T)
cvt <- data.frame(motif_ID=x[,1], s_ref=as.numeric(x[,2]), s_alt=as.numeric(x[,3]))


###
motif0 <- cvt$motif_ID[3]
fn <- "/nfs/rprdata/julong/sc-atac/analyses.2021-02-05/3_motif/2_motif.activities.outs/3_motif.diff.results.rds"
resDiff <- read_rds(fn)%>%mutate(conditions=paste(MCls, contrast, sep="_"))

resDiff%>%filter(gene==motif0, conditions=="Tcell_LPS-DEX")


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
