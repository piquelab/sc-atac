#module load R/test_4.1
###
#setwd("/nfs/rprdata/julong/sc-atac/analyses.2021-02-05/1_seurat")
library(tidyverse)
## library(parallel)
library(data.table)
## library(purrr)
library(GenomicRanges)
library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(Signac)
library(SeuratWrappers)
## library(cicero, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
## library(monocle3)
###
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(viridis)
theme_set(theme_grey())

outdir <- "./1_Merge.outs/"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)

rm(list=ls())

###########################
### 0, generate folders ###
###########################

basefolder <- "/nfs/rprdata/julong/sc-atac/count.SCAIP.2021-01-14/"
expNames <- dir(basefolder,"^SCAIP*")
folders <- paste0(basefolder, expNames, "/", sep="")
ind <- dir.exists(folders)
folders <- folders[ind]
expNames <- expNames[ind]
names(folders) <- expNames


##################################################################
### 1, creating a common peak set from filtered peak.bed files ###
##################################################################

expNames <- names(folders)
peaks <- lapply(expNames, function(ii){
   fn <- paste(folders[ii], "outs/filtered_peak_bc_matrix/peaks.bed", sep="")
   bed <- read.table(fn, col.names=c("chr", "start", "end"))  
   gr <- makeGRangesFromDataFrame(bed)
   gr  
})
x2 <- unlist(GRangesList(peaks))
combined.peaks <- GenomicRanges::reduce(x=x2)
###
peakwidths <- GenomicRanges::width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths<10000&peakwidths>20]
write_rds(combined.peaks, file="./1_Merge.outs/0_combined.peaks.rds")


###
### 3, combined seurat object

### readATAC, quantify peaks and create ATAC object
readATAC <- function(run, peaks){

### load metadata
### extract filtered barcodes from cell ranger
   metafn <- paste(run, "outs/filtered_peak_bc_matrix/barcodes.tsv", sep="")
   meta <- fread(metafn, header=F)
                                                                         
   fragfn <- paste(run, "outs/fragments.tsv.gz", sep="")
### create fragment object
   frags <- CreateFragmentObject(path=fragfn, cells=meta$V1)

### quantify peaks
   counts <- FeatureMatrix(fragments=frags, features=peaks, cells=meta$V1)

### create a seurat object
   assay <- CreateChromatinAssay(counts, fragments=frags)
   atac <- CreateSeuratObject(assay, assay="ATAC")

   atac 
}###


###
### merge multiple objects
combined.peaks <- read_rds("./1_Merge.outs/0_combined.peaks.rds")
expNames <- names(folders)
atac_ls <- future_map(expNames, function(ii){
   cat(ii,"\n")
   atac <- readATAC(run=folders[ii], peaks=combined.peaks)
   atac <- RenameCells(atac, add.cell.id=ii)
   atac
})
###

combined <- merge(atac_ls[[1]], atac_ls[-1], project="sc-atac")

write_rds(combined, file="./1_Merge.outs/1_seurat.merge.rds")


### Add meta data
combined <- read_rds("./1_Merge.outs/1_seurat.merge.rds")

### read meta data across 10 experiments
meta <- map_dfr(expNames, function(ii){
   run <- folders[ii]
   fn <- paste(run, "outs/singlecell.csv", sep="")
   meta <- fread(fn)%>%mutate(barcode=paste(ii, barcode, sep="_"))
   meta
})

x <- combined@meta.data
x <- x%>%mutate(barcode=rownames(x))%>%left_join(meta)
rownames(x) <- x$barcode
combined <- AddMetaData(combined, x)

###add the gene information to the object
annotations <- GetGRangesFromEnsDb(ensdb=EnsDb.Hsapiens.v75)
#seqlevelsStyle(annotations) <- "1000"
genome(annotations) <- "hg19"
Annotation(combined) <- annotations

combined <- combined%>%NucleosomeSignal()%>%TSSEnrichment(fast=F)
combined <- NucleosomeSignal(combined)
#combined <- TSSEnrichment(combined, fast=F)
combined$pct_reads_in_peaks <- combined$peak_region_fragments/combined$passed_filters*100
combined$blacklist_ratio <- combined$blacklist_region_fragments/combined$peak_region_fragments

opfn <- "./1_Merge.outs/1_seurat.merge.rds"
write_rds(combined, file=opfn)            

###
###
 
 
###############
### Summary ###
###############

atac <- read_rds("./1_Merge.outs/1_seurat.merge.rds")
atac$high.tss <- ifelse(atac$TSS.enrichment>2, "High", "Low")

fig0 <- TSSPlot(atac, group.by="high.tss")+
        NoLegend()+
        theme(plot.title=element_text(hjust=0.5))        
figfn <- "./1_Merge.outs/Figure1.tss.png"
png(figfn, width=500, height=400, res=120)
print(fig0)
dev.off() 

##
#atac$nucleosome_group <- ifelse(atac$nucleosome_signal>4, "NS > 4", "NS < 4")
#fig1 <- FragmentHistogram(object=atac, group.by = "nucleosome_group")
#figfn <- "./outs/Figure2.fragment.png"
#png(figfn, width=500, height=400, res=120)
#print(fig1)
#dev.off()

###
fig2 <- VlnPlot(object=atac,
                features=c("pct_reads_in_peaks", "peak_region_fragments", "TSS.enrichment", 
                           "blacklist_ratio", "nucleosome_signal"),
                pt.size=0,
                ncol=3)&
        theme_bw()+
        theme(axis.title=element_blank(),
              axis.text.x=element_blank(),
              legend.position="none",
              plot.title=element_text(hjust=0.5,size=10))
figfn <- "./1_Merge.outs/Figure3.vlnplot.png"
png(figfn, width=700, height=500, res=120)
print(fig2)
dev.off()



