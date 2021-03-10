##
source("../LibraryPackage.R")

outdir <- "./1_ATAC.outs/"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)
rm(list=ls())


#################################
### 1. processing ATAC object ###
################################# 
if (FALSE){

prefix <- "/wsu/home/groups/piquelab/scHOLD/scATAC/counts.scHold.20210301/"

### Read data
### counts
countfn <- paste(prefix, "scHOLD-ATAC/outs/filtered_peak_bc_matrix.h5", sep="") 
counts <- Read10X_h5(countfn)

### objects 
fragfn <- paste(prefix, "scHOLD-ATAC/outs/fragments.tsv.gz", sep="")
chrom_assay <- CreateChromatinAssay(
   counts=counts, 
   sep=c(":", "-"), 
   genome="hg19",
   fragments=fragfn)
   
scHOLD <- CreateSeuratObject(chrom_assay, assay="ATAC")

### Add meta data
metafn <- paste(prefix, "scHOLD-ATAC/outs/singlecell.csv", sep="")
meta2 <- fread(metafn)
###                
x <- scHOLD@meta.data
x <- x%>%mutate(barcode=rownames(x))%>%left_join(meta2)
rownames(x) <- x$barcode
scHOLD <- AddMetaData(scHOLD, x)

###add the gene information to the object
annotations <- GetGRangesFromEnsDb(ensdb=EnsDb.Hsapiens.v75)
genome(annotations) <- "hg19"
Annotation(scHOLD) <- annotations

scHOLD <- scHOLD%>%NucleosomeSignal()%>%TSSEnrichment(fast=F)%>%NucleosomeSignal()

scHOLD$pct_reads_in_peaks <- scHOLD$peak_region_fragments/combined$passed_filters*100
scHOLD$blacklist_ratio <- scHOLD$blacklist_region_fragments/scHOLD$peak_region_fragments
scHOLD$pct_reads_in_peaks <- scHOLD$peak_region_fragments/scHOLD$passed_filters*100
scHOLD$blacklist_ratio <- scHOLD$blacklist_region_fragments/scHOLD$peak_region_fragments
opfn <- "./1_ATAC.outs/1_scHOLD.ATAC.rds"
write_rds(scHOLD, file=opfn)

### clustering analysis
atac <- read_rds("./1_ATAC.outs/1_scHOLD.ATAC.rds")
atac <- atac%>%
        RunTFIDF()%>%
        FindTopFeatures(min.cutoff = "q0")%>%
        RunSVD(n=100)%>%
        RunUMAP(reduction="lsi", dims=2:50, reduction.name="umap.atac", reduction.key="atacUMAP_")
atac2 <- atac%>%
         FindNeighbors(reduction="lsi", dims=2:50)%>%
         FindClusters(resolution=0.3, verbose=FALSE, algorithm=3)
opfn <- "./1_ATAC.outs/2_scHOLD.ATAC.cluster.rds"
write_rds(atac2, file=opfn)

}

if(FALSE){
atac2 <- read_rds("./outs/1.1_scHOLD.ATAC.cluster.rds")         
fig0 <- DimPlot(object=atac2, label=T, raster=F)+NoLegend()+theme_bw()+theme(legend.position="none")
figfn <- "./1_ATAC.outs/Figure2.1_atac.cluster.png"
png(figfn, width=400, height=500, res=120)
print(fig0)
dev.off()  

atac2 <- read_rds("./outs/1.1_scHOLD.ATAC.cluster.rds")
fig0 <- DimPlot(object=atac2, label=F, raster=T)+NoLegend()+
        theme_bw()+
        theme(axis.text=element_blank(), 
              axis.title=element_blank(),
              axis.ticks=element_blank(), legend.position="none")
figfn <- "./1_ATAC.outs/Figure2.2_atac.cluster.png"
png(figfn, width=400, height=500, res=120)
print(fig0)
dev.off()        

###
###
atac2 <- read_rds("./outs/1.1_scHOLD.ATAC.cluster.rds")
atac2$pct_reads_in_peaks <- atac2$peak_region_fragments/atac2$passed_filters*100
atac2$blacklist_ratio <- atac2$blacklist_region_fragments/atac2$peak_region_fragments
fig0 <- VlnPlot(object=atac2,
                features=c("pct_reads_in_peaks", "TSS.enrichment"),
                pt.size=0, group.by="orig.ident",
                ncol=2)&
        theme_bw()+
        theme(axis.title=element_blank(),
              axis.text.x=element_blank(),
              legend.position="none",
              plot.title=element_text(hjust=0.5,size=10))
figfn <- "./1_ATAC.outs/Figure1_vlnplot.png"
png(figfn, width=500, height=250, res=120)
print(fig0)
dev.off()
}
