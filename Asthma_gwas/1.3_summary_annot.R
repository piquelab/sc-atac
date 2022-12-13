##
library(tidyverse)
library(data.table)


### passing argument
args=commandArgs(trailingOnly=T)
if (length(args)>0){
   motifFile <- args[1]
   oneMCl <- args[2]
 }else{
   motifFile <- "splitMotif000"
   oneMCl <- "Bcell"
}
 
outdir <- paste("./torus_input/summary_", oneMCl, "/", sep="")
if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)


###
fn <- paste("./Motif_file/", motifFile, sep="")
motifList <- read.table(fn)$V1


###
summ <- map_dfr(1:length(motifList), function(i){
   ## 
   motif <- motifList[i] 
   cat(i, motif, "\n")    
   ###
   ## fn <- paste(outdir, motif, "_union_torus.annot.gz", sep="") ## for cell-types
   fn <- paste("./torus_input/Motif_", oneMCl, "/",  motif, "_torus.annot.gz", sep="") 
   if ( file.exists(fn)){
      ### 
      annot <- fread(fn, header=T, data.table=F)
      summ2 <- data.table(motif=motif, MCls=oneMCl, num=sum(annot$motif_d==1)) 
      summ2
   }    
})

opfn <- paste(outdir, "summ_", motifFile, ".txt", sep="")
write.table(summ, opfn, quote=F, row.names=F, col.names=F)
### End
