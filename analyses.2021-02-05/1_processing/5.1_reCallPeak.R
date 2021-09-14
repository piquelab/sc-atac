##
library(tidyverse)
## library(parallel)
library(data.table)
## library(purrr)
library(GenomicRanges)
library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(SeuratObject)
library(Signac)
library(SeuratWrappers)
library(cicero, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
library(monocle3)
library(EnsDb.Hsapiens.v75)
###
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(viridis)
theme_set(theme_grey())

###
###
outdir <- "./5.1_reCallPeak.outs/"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)
rm(list=ls())


#####################
### re-call peaks ###
#####################

###
atac <- read_rds("./4.2_Integrate.outs/3_scATAC.annot.rds")
###
peaks <- CallPeaks(atac, group.by="MCls",
   macs2.path="/wsu/home/groups/piquelab/apps/el7/anaconda3python/envs/macs2/bin/macs2")
###
opfn <- "./5.1_reCallPeak.outs/1.1_CallPeak.macs.rds"
write_rds(peaks, opfn)

##
## frag <- Fragments(atac)
## counts <- FeatureMatrix(fragments=frag, features=peaks)
## opfn <- "./5.1_reCallPeak.outs/1.2_counts.rds"
## write_rds(counts, opfn)

####
#### readATAC funtion for read atac data 
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
   assay <- CreateChromatinAssay(counts,
      sep=c(":", "-"), genome="hg19", fragments=frags)
   atac <- CreateSeuratObject(assay, assay="ATAC")
   atac 
}###


##############################
### merge multiple objects ###
##############################

### folders
basefolder <- "/nfs/rprdata/julong/sc-atac/count.SCAIP.2021-01-14/"
expNames <- dir(basefolder,"^SCAIP*")
folders <- paste0(basefolder, expNames, "/", sep="")
ind <- dir.exists(folders)
folders <- folders[ind]
expNames <- expNames[ind]
names(folders) <- expNames

peaks <- read_rds("./5.1_reCallPeak.outs/1.1_CallPeaks.macs.rds")
expNames <- names(folders)
atac_ls <- lapply(expNames, function(ii){
   cat(ii,"\n")
   atac <- readATAC(run=folders[ii], peaks=peaks)
   atac <- RenameCells(atac, add.cell.id=ii)
   atac <- RenameCells(atac, new.names=gsub("-1","",Cells(atac)))
   atac
})
###


combined <- merge(atac_ls[[1]], atac_ls[-1], project="sc-atac")

write_rds(combined, file="./5.1_reCallPeak.outs/2_seurat.merge.rds")


### Add meta data
## combined <- read_rds("./5.1_reCallPeak.outs/2_seurat.merge.rds")
###add the gene information to the object
annotations <- GetGRangesFromEnsDb(ensdb=EnsDb.Hsapiens.v75)
#seqlevelsStyle(annotations) <- "1000"
genome(annotations) <- "hg19"
Annotation(combined) <- annotations

atac <- read_rds("./4.2_Integrate.outs/3_scATAC.annot.rds")
x <- atac@meta.data

combined2 <- subset(combined,cells=Cells(atac))
meta <- combined2@meta.data
metaNew <- cbind(meta,x[,28:37])
combined <- AddMetaData(combined2,metaNew)

opfn <- "./5.1_reCallPeak.outs/3_seurat.annot.rds"
write_rds(combined, opfn)



srun -q primary --mem=300g --time=




    12-23:00:00 --pty bash -l
## combined <- combined%>%NucleosomeSignal()%>%TSSEnrichment(fast=F)
## ## combined <- NucleosomeSignal(combined)
## #combined <- TSSEnrichment(combined, fast=F)
## combined$pct_reads_in_peaks <- combined$peak_region_fragments/combined$passed_filters*100
## combined$blacklist_ratio <- combined$blacklist_region_fragments/combined$peak_region_fragments

## opfn <- "./1_Merge.outs/1_seurat.merge.rds"
## write_rds(combined, file=opfn)  































































