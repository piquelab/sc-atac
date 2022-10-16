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
 
outdir <- paste("./torus_input/Motif_", oneMCl, "/", sep="")
if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)


### Asthma gwas
summ <- fread("./torus_input/Asthma_torus_zval.txt.gz", header=T, data.table=F)

bed <- fread("SCAIP_final_bed.gz", header=F, data.table=F)
bed <- bed%>%filter(V6%in%summ$SNP)

chr_pos <- bed$V6
names(chr_pos) <- as.character(bed$V4)

##
fn <- paste("./Motif_file/", motifFile, sep="")
motifList <- read.table(fn)$V1
###
for (i in 1:length(motifList)){

   motif <- motifList[i] 
   cat(i, motif, "\n")
    
   ###
   fn <- paste("../SNPannotation/4_SNPAnnot.outs/pct_0.02_", oneMCl, "/", motif, "_union_torus.annot.gz", sep="")
   if ( file.exists(fn)){
      ### 
      annot <- fread(fn, header=T, data.table=F) 
      annot2 <- annot%>%filter(SNP%in%bed$V4)%>%mutate(chr_pos_grch38=chr_pos[as.character(SNP)]) 
      annot2 <- annot2%>%dplyr::select(chr_pos_grch38, peaking_d)%>%dplyr::rename(SNP=chr_pos_grch38) 
      ###
      gfn <- gzfile(paste(outdir, motif, "_union_torus.annot.gz", sep=""))
      write.table(annot2, gfn, row.names=F, quote=F, sep="\t", col.names=T)
      ###
   }    
}

### End
