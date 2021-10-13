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
library(ChIPseeker, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
###
library(ggplot2)
library(cowplot, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
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


########################
### annotation peaks ###
########################

###
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


###
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
### compare the closest genes
x1 <- read_rds("./2.2_compareRNAandATAC.outs/1_annot.signac.rds")
x1 <- x1%>%dplyr::rename("peak_region"="query_region")%>%
  dplyr::select(peak_region, gene_id)  
    
x2 <- read_rds("./2.2_compareRNAandATAC.outs/2_annot.ChIPseeker.rds")%>%
   as.data.frame()%>%
   mutate(chr=gsub("chr","",seqnames),
          peak_region=paste(chr,start,end,sep="-"))%>%
  dplyr::select(peak_region, geneId)

x <- x1%>%left_join(x2, by="peak_region")

### compare the genes around regions
x2 <- read_rds("./2.2_compareRNAandATAC.outs/2_annot.ChIPseeker.rds")%>%
   as.data.frame()%>%
   mutate(chr=gsub("chr","",seqnames),
          peak_region=paste(chr,start,end,sep="-"))%>%
   dplyr::select(peak_region, flank_geneIds)
x <- x1%>%left_join(x2, by="peak_region")

geneList <- str_split(x$flank_geneIds, ";")
ii <- map_dbl(1:length(geneList), function(i){
  gene0 <- as.character(x[i, "gene_id"])
  x <- ifelse(gene0%in%geneList[[i]], 1, 0)
  x
})



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
png(figfn, width=480, height=320)
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
### flanking region
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




############################
### if enrich are in DVG ###
############################

### previous identified DEGs
fn <- "/nfs/rprdata/julong/SCAIP/analyses/SCAIP-B1-6_2020.03.23/10_RNA.Variance_output/tmp9/3_phiNew.meta"
resDE <- read.table(fn,header=T)%>%dplyr::filter(qval<0.1, abs(beta)>0.5)
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
      labels=c("0"="Not DVG", "1"="DVG"),
      guide=guide_legend(override.aes=list(size=1)))+
   facet_grid(MCls~contrast, scales="free", labeller=labeller(contrast=lab1))+
   xlab(bquote("Expected"~ -log[10]~"("~italic(plain(P))~")"))+
   ylab(bquote("Observed"~-log[10]~"("~italic(plain(P))~")"))+
   theme_bw()+
   theme(legend.title=element_blank(),strip.text=element_text(size=12))

figfn <- "./2.2_compareRNAandATAC.outs/Figure2.5_qq.DVG.cloest.png"
png(figfn, width=750, height=750, res=120)
print(p1)
dev.off()



##################################
### DEG if is enriched in DARs ###
##################################

### fisher test function
cal.fisher <- function(df){
   ###
 resfisher <- map_dfr(1:nrow(df),function(i){
      dmat <- matrix(as.numeric(df[i,]), 2, 2)
      colnames(dmat) <- c("interest", "not.interest")
      rownames(dmat) <- c("in.DAR", "not.DAR")
      res <- fisher.test(dmat)
      res2 <- data.frame(odds=res$estimate, pval.fisher=res$p.value,
         lower=res$conf.int[1], upper=res$conf.int[2])
      res2
   })
   resfisher
}

### 
df.fisher <- function(res, resDP, peakAnno){
   ### 
   comb <- sort(unique(res$comb))
   df.ls <- lapply(1:length(comb), function(i){
      ii <- comb[i] 
      cell <- gsub("_.*", "", ii)
      contrast <- gsub(".*_", "", ii)
      res2 <- res%>%dplyr::filter(comb==ii) ## dplyr::filter(qval<0.1, abs(beta)>0.5)
      resDP2 <- resDP%>%
         dplyr::filter(comb==ii, p.adjusted<0.1, abs(estimate)>0.5)%>%
         dplyr::left_join(peakAnno, by="peak")
      if( nrow(resDP2)>5){
         sigGene <- res2%>%
            dplyr::filter(qval<0.1, abs(beta)>0.5)%>%
            dplyr::pull(gene)
         ###
         interest.in.DARs <- sum(sigGene%in%resDP2$geneId)
         interest.not.DARs <- length(sigGene)-interest.in.DARs
      ###
         notSig <- setdiff(res2$gene, sigGene)
         not.interest.in.DARs <- sum(notSig%in%resDP2$geneId)
         not.interest.not.DARs <- length(notSig)-not.interest.in.DARs
         df2 <- data.frame("interest.in.DARs"=interest.in.DARs,
            "interest.not.DARs"=interest.not.DARs,
            "not.interest.in.DARs"=not.interest.in.DARs,
            "not.interest.not.DARs"=not.interest.not.DARs)
         df2$cell <- cell
         df2$contrast <- contrast
         df2$comb <- ii
      }else{
        df2 <- NA    
      }
      df2
   })    
   df.ls <- df.ls[!is.na(df.ls)]
   df <- do.call(rbind, df.ls)    
   res <- cal.fisher(df[,1:4])
   res2 <- cbind(df, res)
   ### 
   as.data.frame(res2)
}    
   
    

#################
### read data ###
#################



fn <- "./1.2_DiffPeak.outs/2.0_DESeq.results.rds"
resDP <- read_rds(fn)%>%drop_na(p.value)%>%
   mutate(comb=paste(MCls, contrast, sep="_"))%>%
   dplyr::rename("peak"="gene") 


fn <- "./2.2_compareRNAandATAC.outs/2_annot.ChIPseeker.rds"
peakAnno <- read_rds(fn)%>%as.data.frame()%>%
   mutate(peak=paste(gsub("chr", "", seqnames), start, end, sep="-")) 
peakAnno2 <- peakAnno%>%dplyr::select(peak, geneId, SYMBOL)



### test if DEG is DARs 
fn <- "/nfs/rprdata/julong/SCAIP/analyses/SCAIP-B1-6_2020.03.23/6_DEG.CelltypeNew_output/Filter2/2_meta.rds"
res <- read_rds(fn)%>%drop_na(pval)%>%
   mutate(comb=paste(MCls, contrast, sep="_"))%>%
   as.data.frame()
##
resFisher <- df.fisher(res, resDP, peakAnno2)
opfn <- "./2.2_compareRNAandATAC.outs/3.1_enrich.DEG.csv"
write.csv(resFisher, opfn, row.names=F)


### test if DVGs is DARs
fn <- "/nfs/rprdata/julong/SCAIP/analyses/SCAIP-B1-6_2020.03.23/10_RNA.Variance_output/tmp9/3_phiNew.meta"
res <- read.table(fn, header=T) %>%drop_na(pval)%>%
   mutate(comb=paste(MCls, contrast, sep="_"))

resFisher2 <- df.fisher(res, resDP, peakAnno2)
opfn <- "./2.2_compareRNAandATAC.outs/3.2_enrich.DVG.csv"
write.csv(resFisher2, opfn, row.names=F)




#######################
### error bar plots ###
#######################

col.treat <- c("LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
             "PHA"="#a6cee3", "PHA-DEX"="#1f78b4")
col.MCl <- c("Bcell"="#4daf4a", "Monocyte"="#984ea3",
             "NKcell"="#aa4b56", "Tcell"="#ffaa00")

fn <- "./2.2_compareRNAandATAC.outs/3.1_enrich.DEG.csv"
df <- read.csv(fn, header=T)
df2 <- data.frame(odds=df$odds,
   CI.low=df$lower, CI.high=df$upper,
   comb=gsub("_", ".", df$comb),
   MCls=df$cell, contrast=df$contrast)
##                  
p1 <- ggplot(df2, aes(x=odds, y=comb))+
   geom_errorbarh(aes(xmax=CI.high, xmin=CI.low, colour=MCls),
       size=0.5, height=0.2)+
   geom_point(aes(colour=MCls), shape=19, size=1.5)+
   scale_colour_manual(values=col.MCl)+
   geom_vline(aes(xintercept=1), size=0.25, linetype="dashed")+ 
   xlab("Odds ratio")+
   ggtitle("DEGs are enriched in DARs")+    
   theme_bw()+
   theme(plot.title=element_text(hjust=0.5, size=9),
         axis.title.y=element_blank(),
         axis.title.x=element_text(size=9),
         axis.text=element_text(size=8),
         legend.position="none")

###
fn <- "./2.2_compareRNAandATAC.outs/3.2_enrich.DVG.csv"
df <- read.csv(fn, header=T)
df2 <- data.frame(odds=df$odds,
   CI.low=df$lower, CI.high=df$upper,
   comb=gsub("_", ".", df$comb),
   MCls=df$cell, contrast=df$contrast)
p2 <- ggplot(df2, aes(x=odds, y=comb))+
   geom_errorbarh(aes(xmax=CI.high, xmin=CI.low, colour=MCls),
       size=0.5, height=0.2)+
   geom_point(aes(colour=MCls), shape=19, size=1.5)+
   scale_colour_manual(values=col.MCl)+
   geom_vline(aes(xintercept=1), size=0.25, linetype="dashed")+ 
   xlab("Odds ratio")+
   ggtitle("DVGs are enriched in DARs")+    
   theme_bw()+
   theme(plot.title=element_text(hjust=0.5, size=9),
         axis.title.y=element_blank(),
         axis.title.x=element_text(size=9),
         axis.text=element_text(size=8),
         legend.position="none")

###
figfn <- "./2.2_compareRNAandATAC.outs/Figure3_enrich.DAR.png"
png(figfn, width=680, height=480, res=120)
plot_grid(p1, p2, ncol=2)
dev.off()

