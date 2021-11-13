###
library(tidyverse)
library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(SeuratObject)
library(Signac)
library(SeuratWrappers)
library(SeuratData)

library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Hsapiens.1000genomes.hs37d5, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
library(ChIPseeker, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
###
library(ggplot2)
library(cowplot, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
library(grid)
library(gridExtra)
library(ggExtra)
library(RColorBrewer)
library(ggsci)
library(viridis)
library(ComplexHeatmap)
library(circlize)
theme_set(theme_grey())


outdir <- "./4_motif.outs/"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)

### fisher test
cal.fisher <- function(df){
   ###
   resfisher <- map_dfr(1:nrow(df),function(i){
      dmat <- matrix(as.numeric(df[i,]),2,2)
      colnames(dmat) <- c("interest", "not.interest")
      rownames(dmat) <- c("in.motif", "not.motif")
      res <- fisher.test(dmat, alternative="greater")
      res2 <- data.frame(odds=res$estimate, pval.fisher=res$p.value)
      res2
   })
   resfisher
}


###
enrich.motif <- function(atac, anno, feature, resSig, resDP){
    
   ### filter annotation 
   anno <- anno%>%
      dplyr::filter(!grepl("chrX",seqnames),
                    !grepl("chrY", seqnames))%>%
      dplyr::select(seqnames, start, end, annotation, geneId, SYMBOL, peak)
  ###  
  if ( feature=="Promoter"){ 
     anno <- anno%>%dplyr::filter(grepl("Promoter", annotation))
  }
    
   ### enrich analysis for condition separately
   comb <- sort(unique(resSig$comb))
   enriched.motif <- lapply(1:length(comb), function(i){
##    
      ii <- comb[i]
      cat(i, ii, "\n") 
##
      ii0 <- gsub("_Up|_Down", "", ii) 
      resDP2 <- resDP%>%dplyr::filter(comb==ii0)
      bg.DP <- intersect(resDP2$peak, anno$peak)

### top peaks containing DEG or DVG    
      DEG <- resSig%>%
         dplyr::filter(comb==ii)%>%dplyr::pull(gene)%>%as.character()

      top.DP <- anno%>%
         dplyr::filter(geneId%in%DEG)%>%
         dplyr::pull(peak)
       
      n.in <- length(top.DP) 
      n.not <- length(bg.DP)-n.in

### enrichment analysis       
      if ( length(top.DP)>5){
         enrich2 <- FindMotifs(atac,
            features=top.DP,
            background=bg.DP)
         df <- enrich2%>%dplyr::select(observed, background)%>%
            mutate("interest.in.motif"=observed,
            "interest.not.motif"=n.in-observed,
            "not.interest.in.motif"=background-observed,
            "not.interest.not.motif"=n.not-not.interest.in.motif)
         fisher <- cal.fisher(df[,3:6])
         enrich2 <- cbind(enrich2,fisher)
         enrich2$qvalue.hyper <- p.adjust(enrich2$pvalue, "BH")
         enrich2$qvalue.fisher <- p.adjust(enrich2$pval.fisher, "BH")
         enrich2$comb <- ii
      }else{
         enrich2 <- NA
      }
      enrich2
   })

###
   enriched.motif <- enriched.motif[!is.na(enriched.motif)]
   enriched.motif <- do.call(rbind, enriched.motif)
   enriched.motif
} ###    


#################
### read data ###
#################

###
### read atac data
fn <- "./2_motif.activities.outs/1_scATAC.motifActivities.rds"
atac <- read_rds(fn)

###
### peak annotation for promoter
fn <- "../2_Differential/2.2_compareRNAandATAC.outs/2_annot.ChIPseeker.rds"
anno <- read_rds(fn)%>%as.data.frame()%>%
   mutate(peak=paste(gsub("chr","",seqnames), start, end, sep="-")) 


###
### results of differential peak
fn <- "../2_Differential/1.2_DiffPeak.outs/2.0_DESeq.results.rds"
resDP <- read_rds(fn)%>%as.data.frame()%>%
   dplyr::filter(MCls!="DC")%>%
   dplyr::rename("peak"="gene")%>%
   mutate(comb=paste(MCls, contrast, sep="_"))


########################
### motif enrich DEG ###
########################

### results of DEG
fn <- "/nfs/rprdata/julong/SCAIP/analyses/SCAIP-B1-6_2020.03.23/6_DEG.CelltypeNew_output/Filter2/2_meta.rds"
resSig <- read_rds(fn)%>%as.data.frame()%>%
   mutate(comb=paste(MCls, contrast, sep="_"))%>% 
   dplyr::filter(qval<0.1, abs(beta)>0.5) 
enrich <- enrich.motif(atac, anno2, resSig, resDP)
write_rds(enrich, "./4_motif.outs/1.1_motif.DEG.rds")


###############
### dotplot ###
###############

##label
clust <- paste(rep(c("LPS", "PHA", "LPS-DEX", "PHA-DEX"), each=4),
   rep(c("Bcell", "Monocyte", "NKcell", "Tcell"), times=4), sep=".")
###
clust2 <- paste(rep(c("A","B","C","D"),each=4), rep(1:4,times=4), sep="")
names(clust2) <- clust
###
lab2 <- gsub("-", "+", clust)
names(lab2) <- clust2


enrich <- read_rds("./4_motif.outs/1.1_motif.DEG.rds")%>%
   mutate(MCls=gsub("_.*", "", comb),
          contrast=gsub(".*_", "", comb),
          comb2=paste(contrast, MCls, sep="."),
          newCluster=clust2[as.character(comb2)])

enrich2 <- enrich%>%
   dplyr::filter(fold.enrichment>1, qvalue.hyper<0.1)%>%
   group_by(newCluster)%>%
   top_n(n=5, wt=fold.enrichment)%>%ungroup()

topmotif <- enrich2%>%dplyr::pull(motif)

enrich3 <- enrich%>%
   dplyr::filter(motif%in%topmotif, fold.enrichment>1, qvalue.hyper<0.1) 

###
p <- ggplot(enrich3, aes(x=newCluster, y=motif.name))+
   geom_point(aes(size=fold.enrichment, colour=qvalue.hyper))+
   scale_colour_gradient(name="FDR",
      low="blue", high="red", na.value=NA, trans="reverse", n.breaks=5,
      guide=guide_colourbar(order=1))+
   scale_size_binned("fold enrichment",
      guide=guide_bins(show.limits=T, axis=T,
      axis.show=arrow(length=unit(1.5,"mm"), ends="both"), order=2),
      n.breaks=4)+
   scale_x_discrete(labels=lab2)+
   theme_bw()+
   theme(axis.title=element_blank(),
         axis.text.x=element_text(angle=60, hjust=1, size=10),
         axis.text.y=element_text(size=8))
###
figfn <- "./4_motif.outs/Figure1.1_dotplot.png"
png(figfn, width=1000, height=1000, res=120)
print(p)
dev.off()


##############################################################
###  motifs are enriched in DEG for up and down separately ###
##############################################################

##

### results of DEG
fn <- "/nfs/rprdata/julong/SCAIP/analyses/SCAIP-B1-6_2020.03.23/6_DEG.CelltypeNew_output/Filter2/2_meta.rds"
res <- read_rds(fn)
resSig <- res%>%as.data.frame()%>%
   mutate(direction=ifelse(beta>0, "Up", "Down"),
          comb=paste(MCls, contrast, direction, sep="_"))%>% 
   dplyr::filter(qval<0.1, abs(beta)>0.5)
##
genetest <- unique(res$gene)
geneList <- str_split(anno$flank_geneIds, ";")
is_gene <- map_lgl(geneList, ~any((.x)%in%genetest) )
anno$is_gene <- is_gene 
anno2 <- anno%>%dplyr::filter(is_gene)

feature0 <- "Allpeaks"
i <- 3
enrich <- enrich.motif(atac, anno2, feature=feature0, resSig, resDP)
opfn <- paste("./4_motif.outs/1.", i, "_motif.DEG.direction.", feature0, ".rds", sep="")
write_rds(enrich, opfn)



#################
### dot plots ###
#################

## re-label
clust <- paste(rep(c("LPS", "PHA", "LPS-DEX", "PHA-DEX"), each=4),
                 rep(c("Bcell", "Monocyte", "NKcell", "Tcell"), times=4), sep=".")
clust <- paste(rep(clust,times=2), rep(c("Up","Down"), each=16), sep=".")
###
ii <- paste(rep(c("A","B","C","D"),each=4),rep(1:4,times=4), sep="")
clust2 <- paste(rep(c("X","Y"),each=16), rep(ii,times=2), sep=".")
names(clust2) <- clust

lab2 <- gsub("-", "+", clust)
names(lab2) <- clust2

## read data for dot plots
feature0 <- "Allpeaks"
i <- 3
fn <- paste("./4_motif.outs/1.", i, "_motif.DEG.direction.", feature0, ".rds", sep="")
enrich <- read_rds(fn)%>%
   mutate(MCls=gsub("_.*", "", comb),
          contrast=gsub(".*[le]_|_[UD].*", "", comb),
          direction=gsub(".*_", "", comb),
          comb2=paste(contrast, MCls, direction, sep="."),
          newCluster=clust2[as.character(comb2)])%>%as.data.frame()

## top 5 motifs
enrich2 <- enrich%>%
   dplyr::filter(fold.enrichment>1, qvalue.hyper<0.1)%>%
   group_by(newCluster)%>%
   top_n(n=5, wt=fold.enrichment)%>%ungroup()

topmotif <- enrich2%>%dplyr::pull(motif)
topmotif <- unique(topmotif)

enrich3 <- enrich%>%
   dplyr::filter(motif%in%topmotif, fold.enrichment>1, qvalue.hyper<0.1) 

###
p <- ggplot(enrich3, aes(x=newCluster, y=motif.name))+
   geom_point(aes(size=fold.enrichment, colour=qvalue.hyper))+
   scale_colour_gradient(name="FDR",
      low="blue", high="red", na.value=NA, trans="reverse", n.breaks=5,
      guide=guide_colourbar(order=1))+
   scale_size_binned("fold enrichment",
      guide=guide_bins(show.limits=T, axis=T,
      axis.show=arrow(length=unit(1.5,"mm"), ends="both"), order=2),
      n.breaks=4)+
   scale_x_discrete(labels=lab2)+
   theme_bw()+
   theme(axis.title=element_blank(),
         axis.text.x=element_text(angle=60, hjust=1, size=10),
         axis.text.y=element_text(size=8))
###
figfn <- paste("./4_motif.outs/Figure1.", i, "_dotplot.direction.", feature0, ".png", sep="")
png(figfn, width=1000, height=1000, res=120)
print(p)
dev.off()



###################################################################
###  motifs are enriched in DEG for up and down separately, NB  ###
###################################################################

### results of DEG
fn <- "/nfs/rprdata/julong/SCAIP/analyses/SCAIP-B1-6_2020.03.23/10_RNA.Variance_output/tmp9/4_mu.meta"
res <- read.table(fn, header=T)
resSig <- res%>%as.data.frame()%>%
   mutate(direction=ifelse(beta>0, "Up", "Down"),
          comb=paste(MCls, contrast, direction, sep="_"))%>% 
   dplyr::filter(qval<0.1, abs(beta)>0.5)
##
genetest <- unique(res$gene)
geneList <- str_split(anno$flank_geneIds, ";")
is_gene <- map_lgl(geneList, ~any((.x)%in%genetest) )
anno$is_gene <- is_gene 
anno2 <- anno%>%dplyr::filter(is_gene)


feature0 <- "Promoter"
i <- 2
enrich <- enrich.motif(atac, anno2, feature=feature0, resSig, resDP)
opfn <- paste("./4_motif.outs/2.", i, "_motif.DEG.direction.", feature0, ".rds", sep="")
write_rds(enrich, opfn)



#################
### dot plots ###
#################
### labels
clust <- paste(rep(c("LPS", "PHA", "LPS-DEX", "PHA-DEX"), each=4),
                 rep(c("Bcell", "Monocyte", "NKcell", "Tcell"), times=4), sep=".")
clust <- paste(rep(clust,times=2), rep(c("Up","Down"), each=16), sep=".")
###
ii <- paste(rep(c("A","B","C","D"),each=4),rep(1:4,times=4), sep="")
clust2 <- paste(rep(c("X","Y"),each=16), rep(ii,times=2), sep=".")
names(clust2) <- clust

lab2 <- gsub("-", "+", clust)
names(lab2) <- clust2

### read data for plot
feature0 <- "Allpeaks"
i <- 3
fn <- paste("./4_motif.outs/2.", i, "_motif.DEG.direction.", feature0, ".rds", sep="")
enrich <- read_rds(fn)%>%
   mutate(MCls=gsub("_.*", "", comb),
          contrast=gsub(".*[le]_|_[UD].*", "", comb),
          direction=gsub(".*_", "", comb),
          comb2=paste(contrast, MCls, direction, sep="."),
          newCluster=clust2[as.character(comb2)])%>%as.data.frame()

### top motifs
enrich2 <- enrich%>%
   ## dplyr::filter(fold.enrichment>1, qvalue.hyper<0.1)%>%
   group_by(newCluster)%>%
   top_n(n=5, wt=fold.enrichment)%>%ungroup()

topmotif <- enrich2%>%dplyr::pull(motif)
topmotif <- unique(topmotif)

enrich3 <- enrich%>%
   dplyr::filter(motif%in%topmotif, fold.enrichment>1, qvalue.hyper<0.1) 

###
p <- ggplot(enrich3, aes(x=newCluster, y=motif.name))+
   geom_point(aes(size=fold.enrichment, colour=qvalue.hyper))+
   scale_colour_gradient(name="FDR",
      low="blue", high="red", na.value=NA, trans="reverse", n.breaks=5,
      guide=guide_colourbar(order=1))+
   scale_size_binned("fold enrichment",
      guide=guide_bins(show.limits=T, axis=T,
      axis.show=arrow(length=unit(1.5,"mm"), ends="both"), order=2),
      n.breaks=4)+
   scale_x_discrete(labels=lab2)+
   theme_bw()+
   theme(axis.title=element_blank(),
         axis.text.x=element_text(angle=60, hjust=1, size=10),
         axis.text.y=element_text(size=8))
###
figfn <- paste("./4_motif.outs/Figure2.", i, "_dotplot.DEG.direction.", feature0, ".png", sep="")
png(figfn, width=1000, height=1000, res=120)
print(p)
dev.off()




##############################################################
###  motifs are enriched in DVG for up and down separately ###
##############################################################

### results of DVG
fn <- "/nfs/rprdata/julong/SCAIP/analyses/SCAIP-B1-6_2020.03.23/10_RNA.Variance_output/tmp9/3_phiNew.meta"
res <- read.table(fn, header=T)
resSig <- res%>%as.data.frame()%>%
   mutate(direction=ifelse(beta>0, "Up", "Down"),
          comb=paste(MCls, contrast, direction, sep="_"))%>% 
   dplyr::filter(qval<0.1, abs(beta)>0.5)

##
genetest <- unique(res$gene)
geneList <- str_split(anno$flank_geneIds, ";")
is_gene <- map_lgl(geneList, ~any((.x)%in%genetest) )
anno$is_gene <- is_gene 
anno2 <- anno%>%dplyr::filter(is_gene)

feature0 <- "Allpeaks"
i <- 3
enrich <- enrich.motif(atac, anno2, feature=feature0, resSig, resDP)

opfn <- paste("./4_motif.outs/3.", i, "_motif.DVG.direction.", feature0, ".rds", sep="")
write_rds(enrich, opfn)


#################################
### dot plot for enrich motif ###
#################################

### re-labels
clust <- paste(rep(c("LPS", "PHA", "LPS-DEX", "PHA-DEX"), each=4),
                 rep(c("Bcell", "Monocyte", "NKcell", "Tcell"), times=4), sep=".")
clust <- paste(rep(clust,times=2), rep(c("Up","Down"), each=16), sep=".")
###
ii <- paste(rep(c("A","B","C","D"),each=4),rep(1:4,times=4), sep="")
clust2 <- paste(rep(c("X","Y"),each=16), rep(ii,times=2), sep=".")
names(clust2) <- clust

lab2 <- gsub("-", "+", clust)
names(lab2) <- clust2

### read data for dot plots
feature0 <- "Allpeaks"
i <- 3
fn <- paste("./4_motif.outs/3.", i, "_motif.DVG.direction.", feature0, ".rds", sep="")
enrich <- read_rds(fn)%>%
   mutate(MCls=gsub("_.*", "", comb),
          contrast=gsub(".*[le]_|_[UD].*", "", comb),
          direction=gsub(".*_", "", comb),
          comb2=paste(contrast, MCls, direction, sep="."),
          newCluster=clust2[as.character(comb2)])%>%as.data.frame()

### top 10 motifs
enrich2 <- enrich%>%
   ## dplyr::filter(fold.enrichment>1, qvalue.hyper<0.1)%>%
   group_by(newCluster)%>%
   top_n(n=5, wt=fold.enrichment)%>%ungroup()

topmotif <- enrich2%>%dplyr::pull(motif)
topmotif <- unique(topmotif)

enrich3 <- enrich%>%
   dplyr::filter(motif%in%topmotif, fold.enrichment>1, qvalue.hyper<0.1) 

###
p <- ggplot(enrich3, aes(x=newCluster, y=motif.name))+
   geom_point(aes(size=fold.enrichment, colour=qvalue.hyper))+
   scale_colour_gradient(name="FDR",
      low="blue", high="red", na.value=NA, trans="reverse", n.breaks=5,
      guide=guide_colourbar(order=1))+
   scale_size_binned("fold enrichment",
      guide=guide_bins(show.limits=T, axis=T,
      axis.show=arrow(length=unit(1.5,"mm"), ends="both"), order=2),
      n.breaks=4)+
   scale_x_discrete(labels=lab2)+
   theme_bw()+
   theme(axis.title=element_blank(),
         axis.text.x=element_text(angle=60, hjust=1, size=10),
         axis.text.y=element_text(size=8))
###
figfn <- paste("./4_motif.outs/Figure3.", i, "_dotplot.DVG.direction.", feature0, ".png", sep="")
png(figfn, width=1000, height=1000, res=120)
print(p)
dev.off()
