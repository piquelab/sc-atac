###
library(Matrix)
library(tidyverse)
library(qqman)
library(qvalue)
##
library(DESeq2)
library(biobroom)
library(ashr)
library(GenomicRanges)
library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(Signac)
### annotation required package
library(clusterProfiler)
library(org.Hs.eg.db)
library(annotables)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(EnsDb.Hsapiens.v75)
library(ChIPseeker,lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
###
library(ggplot2)
library(cowplot,lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
library(grid)
library(gridExtra)
library(ggExtra)
library(RColorBrewer)
library(ggsci)
library(RColorBrewer)
library(viridis)
theme_set(theme_grey())


outdir <- "./2.2_compareRNAandATAC.outs/"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F)


### annotation database
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
edb <- EnsDb.Hsapiens.v75
seqlevelsStyle(edb) <- "UCSC"

atac <- read_rds("../1_processing/5.1_reCallPeak.outs/3_scATAC.annot.rds")
anno <- Annotation(atac)

###
### peak annotation from signac 
peakAnno <- ClosestFeature(atac, regions=granges(atac))
opfn <- "./2.2_compareRNAandATAC.outs/1_annot.signac.rds"
write_rds(peakAnno, opfn)
## peakAnno <- peakAnno%>%dplyr::select(gene_name,query_region)


### peak annotation from ChIPseeker
x <- as.data.frame(granges(atac))%>%
   mutate(seqnames=paste("chr",seqnames,sep=""))
peak <- makeGRangesFromDataFrame(x)
peakAnno <- annotatePeak(peak, tssRegion=c(-3000,3000), TxDb=edb,
   addFlankGeneInfo=T, flankDistance=100000, annoDb="org.Hs.eg.db")
###
opfn <- "./2.2_compareRNAandATAC.outs/2_annot.ChIPseeker.rds"
write_rds(peakAnno, opfn)


###
### compare signac and ChIPseeker 
x1 <- read_rds("./2.2_compareRNAandATAC.outs/1_annot.signac.rds")
x1 <- x1%>%dplyr::rename("peak_region"="query_region")%>%
  dplyr::select(peak_region, gene_id)  
    
x2 <- read_rds("./2.2_compareRNAandATAC.outs/2_annot.ChIPseeker.rds")%>%
   as.data.frame(peakAnno)%>%
   mutate(chr=gsub("chr","",seqnames),
          peak_region=paste(chr,start,end,sep="-"))%>%
  dplyr::select(peak_region,geneId)

x <- x1%>%left_join(x2,by="peak_region")


#################################################
### if enriched in DEG from signac annotation ###
#################################################

res <- read_rds("./1.2_DiffPeak.outs/2.0_DESeq.results.rds")%>%as.data.frame()
res <- res%>%dplyr::rename("peak_region"="gene")%>%drop_na(p.value)

peakAnno <- read_rds("./2.2_compareRNAandATAC.outs/1_annot.signac.rds")
peakAnno <- peakAnno%>%
   dplyr::rename("peak_region"="query_region")%>%
   dplyr::select(gene_id, peak_region)

res <- res%>%
   left_join(peakAnno, by="peak_region")%>%
   mutate(comb=paste(MCls, contrast, sep="_"))


### previous identified DEGs
fn <- "/nfs/rprdata/julong/SCAIP/analyses/SCAIP-B1-6_2020.03.23/6_DEG.CelltypeNew_output/Filter2/2_meta.rds"
resDE <- read_rds(fn)%>%dplyr::filter(qval<0.1, abs(beta)>0.5)%>%
   mutate(comb=paste(MCls, contrast, sep="_"))


comb <- sort(unique(resDE$comb))
dfNew <- map_dfr(comb, function(ii){
  res2 <- res%>%dplyr::filter(comb==ii)
  DEG <- resDE%>%dplyr::filter(comb==ii)%>%dplyr::pull(gene)
  res2 <- res2%>%mutate(is_DEG=ifelse(gene_id%in%DEG,1,0))
  ###
  dx <- map_dfr(c(0,1),function(i){
     di <- res2%>%dplyr::filter(is_DEG==i)%>%arrange(p.value)
     ngene <- nrow(di)
     di <- di%>%mutate(observed=-log10(p.value), expected=-log10(ppoints(ngene)))
     di
  })
  dx
})


###
###
lab1 <- c("LPS"="LPS", "LPS-DEX"="LPS+DEX", "PHA"="PHA", "PHA-DEX"="PHA+DEX")
p1 <- ggplot(dfNew, aes(x=expected,y=observed, colour=factor(is_DEG)))+
   ggrastr::rasterise(geom_point(size=0.3),dpi=300)+
   geom_abline(colour="red")+
   scale_colour_manual(values=c("0"="grey40","1"="green"),
      labels=c("0"="Not DEG", "1"="DEG"),
      guide=guide_legend(override.aes=list(size=1)))+
   facet_grid(MCls~contrast, scales="free", labeller=labeller(contrast=lab1))+
   xlab(bquote("Expected"~ -log[10]~"("~italic(plain(P))~")"))+
   ylab(bquote("Observed"~-log[10]~"("~italic(plain(P))~")"))+
   theme_bw()+
   theme(legend.title=element_blank(),strip.text=element_text(size=12))

figfn <- "./2.2_compareRNAandATAC.outs/Figure1.1_qq.signac.png"
png(figfn, width=750, height=750, res=120)
print(p1)
dev.off()



#######################################
### peak annotation from ChIPseeker ###
#######################################

peakAnno <- read_rds("./2.2_compareRNAandATAC.outs/2_annot.ChIPseeker.rds")

##
figfn <- "./2.2_compareRNAandATAC.outs/Figure2.1_annot.pie.png"
png(figfn)
plotAnnoPie(peakAnno)
dev.off()

###
figfn <- "./2.2_compareRNAandATAC.outs/Figure2.2_annot.ven.png"
png(figfn, width=700, height=500)
vennpie(peakAnno)
dev.off()

###
peakAnno <- as.data.frame(peakAnno)
peakAnno <- peakAnno%>%
   mutate(chr=gsub("chr","",seqnames),
          peak_region=paste(chr,start,end,sep="-"))%>%
   dplyr::select(peak_region, geneId, flank_geneIds) 


###
res <- read_rds("./1.2_DiffPeak.outs/2.0_DESeq.results.rds")%>%as.data.frame()
res <- res%>%dplyr::rename("peak_region"="gene")%>%drop_na(p.value)

res <- res%>%
   left_join(peakAnno, by="peak_region")%>%
   mutate(comb=paste(MCls, contrast, sep="_"))


### previous identified DEGs
fn <- "/nfs/rprdata/julong/SCAIP/analyses/SCAIP-B1-6_2020.03.23/6_DEG.CelltypeNew_output/Filter2/2_meta.rds"
resDE <- read_rds(fn)%>%dplyr::filter(qval<0.1, abs(beta)>0.5)
resDE <- resDE%>%
   mutate(comb=paste(MCls, contrast, sep="_"))


###
### qq for cloest genes
comb <- sort(unique(resDE$comb))
dfNew <- map_dfr(comb, function(ii){
  res2 <- res%>%dplyr::filter(comb==ii)
  DEG <- resDE%>%dplyr::filter(comb==ii)%>%dplyr::pull(gene)
  res2 <- res2%>%mutate(is_DEG=ifelse(geneId%in%DEG,1,0))  
  ###
  dx <- map_dfr(c(0,1),function(i){
     di <- res2%>%dplyr::filter(is_DEG==i)%>%arrange(p.value)
     ngene <- nrow(di)
     di <- di%>%mutate(observed=-log10(p.value), expected=-log10(ppoints(ngene)))
     di
  })
  dx
})

###
###
lab1 <- c("LPS"="LPS", "LPS-DEX"="LPS+DEX", "PHA"="PHA", "PHA-DEX"="PHA+DEX")
p1 <- ggplot(dfNew, aes(x=expected,y=observed, colour=factor(is_DEG)))+
   ggrastr::rasterise(geom_point(size=0.3),dpi=300)+
   geom_abline(colour="red")+
   scale_colour_manual(values=c("0"="grey40","1"="green"),
      labels=c("0"="Not DEG", "1"="DEG"),
      guide=guide_legend(override.aes=list(size=1)))+
   facet_grid(MCls~contrast, scales="free", labeller=labeller(contrast=lab1))+
   xlab(bquote("Expected"~ -log[10]~"("~italic(plain(P))~")"))+
   ylab(bquote("Observed"~-log[10]~"("~italic(plain(P))~")"))+
   theme_bw()+
   theme(legend.title=element_blank(),strip.text=element_text(size=12))

figfn <- "./2.2_compareRNAandATAC.outs/Figure2.3_qq.DEG.cloest.png"
png(figfn, width=750, height=750, res=120)
print(p1)
dev.off()



###
### 
### generate data for qq plot
comb <- sort(unique(resDE$comb))
dfNew <- map_dfr(comb, function(ii){
  res2 <- res%>%dplyr::filter(comb==ii)
  DEG <- resDE%>%dplyr::filter(comb==ii)%>%dplyr::pull(gene)
  geneList <- str_split(res2$flank_geneIds, ";")
  is_DEG <- map_dbl(geneList, function(i){
     x <- ifelse(any(i%in%DEG),1,0)
     x
  })
  res2$is_DEG <- is_DEG
  
  ###
  dx <- map_dfr(c(0,1),function(i){
     di <- res2%>%dplyr::filter(is_DEG==i)%>%arrange(p.value)
     ngene <- nrow(di)
     di <- di%>%mutate(observed=-log10(p.value), expected=-log10(ppoints(ngene)))
     di
  })
  dx
})


###
lab1 <- c("LPS"="LPS", "LPS-DEX"="LPS+DEX", "PHA"="PHA", "PHA-DEX"="PHA+DEX")
p2 <- ggplot(dfNew, aes(x=expected,y=observed, colour=factor(is_DEG)))+
   ggrastr::rasterise(geom_point(size=0.3),dpi=300)+
   geom_abline(colour="red")+
   scale_colour_manual(values=c("0"="grey40","1"="green"),
      labels=c("0"="Not DEG", "1"="DEG"),
      guide=guide_legend(override.aes=list(size=1)))+
   facet_grid(MCls~contrast, scales="free", labeller=labeller(contrast=lab1))+
   xlab(bquote("Expected"~ -log[10]~"("~italic(plain(P))~")"))+
   ylab(bquote("Observed"~-log[10]~"("~italic(plain(P))~")"))+
   theme_bw()+
   theme(legend.title=element_blank(),strip.text=element_text(size=12))

figfn <- "./2.2_compareRNAandATAC.outs/Figure2.4_qq.DEG.png"
png(figfn, width=750, height=750, res=120)
print(p2)
dev.off()


