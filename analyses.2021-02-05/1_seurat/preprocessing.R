#module load R/test_4.0.3
###
setwd("./1_seurat/")
source("../LibraryPackage.R")

outdir <- "./outs/"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)

################################################
### 1, generate folders containing h5ad data ###
################################################

basefolder <- "/nfs/rprdata/julong/sc-atac/count.SCAIP.2021-01-14/"
expNames <- dir(basefolder,"^SCAIP*")
folders <- paste0(basefolder, expNames, "/", sep="")
ind <- dir.exists(folders)
folders <- folders[ind]
expNames <- expNames[ind]
names(folders) <- expNames

###
###
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
write_rds(combined.peaks, file="./outs/combined.peaks.rds")

### read ATAC project
readATAC <- function(run, peaks, threshold){

metafn <- paste(run, "outs/singlecell.csv", sep="")
meta <- fread(metafn)%>%filter(passed_filters>threshold)

fragfn <- paste(run, "outs/fragments.tsv.gz", sep="")
frags <- CreateFragmentObject(path=fragfn, cells=meta$barcode)

### quantify peaks
counts <- FeatureMatrix(fragments=frags, features=combined.peaks, cells=meta$barcode) 
}


