###
library(tidyverse)
library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(SeuratObject)
library(Signac)
library(SeuratWrappers)
library(SeuratData)
library(qqman)
library(GenomicRanges)

## library(JASPAR2020)
## library(TFBSTools)
## library(BSgenome.Hsapiens.UCSC.hg19)
## library(BSgenome.Hsapiens.1000genomes.hs37d5, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
## library(ChIPseeker, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
## ###
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


###
#### cell-type active peaks and related motif



###############################
### cell-type active peaks  ###
###############################

rm(list=ls())

### motif

fn <- "./1_motif.outs/1_scATAC.motif.rds" 
atac <- read_rds(fn)
motif <- Motifs(atac)
###
motif.name <- unlist(motif@motif.names)

X <- motif@data
## ##
## opfn <- "./1.3_motif.outs/0_motif.rds"
## write_rds(X, opfn)

x <- granges(atac)
xrange <- ranges(x)
count <- atac@assays$ATAC@counts
anno <- data.frame(rn=rownames(count), rnz=rowSums(count), chr=as.character(seqnames(x)))
autosome <- as.character(1:22)
annoSel <- anno%>%dplyr::filter(rnz>20)
peakAll <- annoSel$rn
                   
## pfm <- GetMotifData(object=motif, slot="pwm")
## motif <- SetMotifData(object=motif, slot="pwm", new.data=pfm)
## differential peaks

###
###

## grid values
pct_grid <- c(0.1, 0.05, 0.02, 0.01)
for(k in 2:length(pct_grid)){
###
pct0 <- pct_grid[k]    
outdir <- paste("./1.3_motif.outs/pct_", pct0, "/", sep="")
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)
##
cat("pct", pct0, sep="")    

    
MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
###
dz <- map_dfc(1:length(MCls), function(i){
###
   oneMCl <- MCls[i]
   atac2 <- subset(atac, subset=MCls==oneMCl)
   count2 <- atac2@assays$ATAC@counts
   rpz <- rowMeans(count2>0)
###
   peakSel <- rownames(count2)[rpz>pct0]
   z <- ifelse(peakAll%in%peakSel, 1, 0)
   z
})
dz <- cbind(peakAll, dz)
names(dz)[2:5] <- paste("d_", MCls, sep="")
rownames(dz) <- dz$peakAll
###
opfn <- paste(outdir, "1_cell-type_active.peaks.rds", sep="")
write_rds(dz, opfn)
}  ###  



#########################################################
### fit glmnet model to select cell-type active motif ###
#########################################################


rm(list=ls())

library(glmnet)
library(Matrix)
library(doMC)
library(parallel)
library(foreach)
registerDoMC(cores=10)


###
### read motif data 
fn <-  "./1.3_motif.outs/0_motif.rds"
motif <- read_rds(fn)
MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")


### loop by a grid of pct
pct_grid <- c(0.1, 0.05, 0.02, 0.01)
for (i in 2:length(pct_grid)){
###
pct0 <- pct_grid[i]

cat("pct", pct0, "\n")
    
outdir <- paste("./1.3_motif.outs/pct_", pct0, "/", sep="")

fn <- paste(outdir, "1_cell-type_active.peaks.rds", sep="")
dz <- read_rds(fn)

###
### loop by cell-type
## B <- lapply(1:4, function(i){
## ###
   k <- 4
   oneMCl <- MCls[k]
###    
   Y <- dz[,k+1]
   X <- motif[dz$peakAll,]
## fit lasso
   system.time(cvfit <- cv.glmnet(X, Y, family="binomial", type.measure="class"))
#
   opfn <- paste(outdir, "2_", k, "_", oneMCl, ".glmnet.rds", sep="")
   write_rds(cvfit, opfn)

##    
   b <- coef(cvfit, s="lambda.1se")
   b <- as.matrix(b)
## output
   opfn <- paste(outdir, "2_", k, "_", oneMCl, ".coef.rds", sep="")
   write_rds(b, opfn)
} ###



## })
## B <- do.call(cbind, B)

###
## colnames(B) <- paste("coef_", MCls, sep="")
##
## opfn <- paste(outdir, "2_coef_motif.rds", sep="")
## write_rds(B, opfn)



#############
### plots ###
#############

##
## cvplot
pct0 <- 0.02
outdir <- paste("./1.3_motif.outs/pct_", pct0, "/", sep="")
oneMCl <- "Tcell"
##
fn <- paste(outdir, "2_4_Tcell.glmnet.rds", sep="")
cvfit <- read_rds(fn)


### a grid of lambda plots
figfn <- paste(outdir, "Figure0.4_Tcell.lambda.png", sep="")
png(figfn, width=520, height=380, res=100)
plot(cvfit)
dev.off()



###
### coef plot
###
b <- as.matrix(coef(cvfit,s="lambda.min"))
b2 <- b[-1,]
df <- data.frame(x=1:length(b2), y=b2)
p1 <- ggplot(df, aes(x,y))+
   geom_point(color="red", size=1)+
   geom_segment(aes(xend=x, yend=0), color="blue", size=0.5)+
   xlab("motif")+ylab("coefficients")+ 
   ggtitle(paste(oneMCl, "(lambda.min)", sep=""))+ 
   theme_bw()+
   theme(plot.title=element_text(hjust=0.5))

###
b <- as.matrix(coef(cvfit,s="lambda.1se"))
b2 <- b[-1,]
df <- data.frame(x=1:length(b2), y=b2)
p2 <- ggplot(df, aes(x,y))+
   geom_point(color="red", size=1)+
   geom_segment(aes(xend=x, yend=0), color="blue", size=0.5)+
   xlab("motif")+ylab("coefficients")+ 
   ggtitle(paste(oneMCl, "(lambda.1se)", sep=""))+ 
   theme_bw()+
   theme(plot.title=element_text(hjust=0.5))   

##
###
figfn <- paste(outdir, "Figure1.4", "_", oneMCl, ".coef.png", sep="")
png(figfn, width=480, height=550, res=120)
plot_grid(p1, p2, nrow=2)
dev.off()



#####################
### select motifs ###
#####################

i <- 4
oneMCl <- "Tcell"

pct_grid <- c(0.1, 0.05, 0.02, 0.01)
for ( k in 2:length(pct_grid)){

    ###
    pct0 <- pct_grid[k]
    dirfn <- paste("./1.3_motif.outs/pct_", pct0, "/", sep="")
    fn <- paste(dirfn, "2_4_Tcell.coef.rds", sep="")
    b <- read_rds(fn)

### se
   bb <- map_dfc(1:20, function(i){
      ##
      fn <- paste(dirfn, "4_Tcell_Bootstrap_", i, "_coef_motif.rds", sep="")
      b <- read_rds(fn) 
      b
   })
###
   bse <- sqrt(apply(bb, 1, var))

  ###  
  res <- data.frame(motif.id=as.character(rownames(b)), coef=b[,1], "bse"=bse)%>%
      mutate(z=coef/bse, pval=pnorm(z, lower.tail=F), "motif.name"=motif.name[motif.id])
  opfn <- paste(dirfn, "2_4_Tcell.zscore.txt", sep="")
  write.table(res, opfn, row.names=F, col.names=T, quote=F, sep="\t")  
  ##
  ## res2 <- res%>%filter(coef>0, pval<0.05)
  ## ### output  
  ## opfn <- paste(dirfn, "2_4_Tcell.selectMotif.txt", sep="")
  ## write.table(res, opfn, row.names=F, col.names=T, quote=F, sep="\t")
    
}

##
## num <- sapply(2:length(pct_grid), function(i){
##     pct0 <- pct_grid[i]
##     fn <- paste("1.3_motif.outs/pct_", pct0, "/2_4_Tcell.selectMotif.txt", sep="")
##     x <- read.table(fn,header=T)
##     nrow(x)
## })


###
### z score plots
## pct0 <- 0.02
## dirfn <- paste("./1.3_motif.outs/pct_", pct0, "/", sep="")
## fn <- paste(dirfn, "2_4_Tcell.zscore.txt", sep="")
## df <- read.table(fn, header=T)
## df <- df[-1,]

## df2 <- data.frame(x=1:nrow(df), y=-log10(df$pval))
## p2 <- ggplot(df2, aes(x,y))+
##    geom_point(color="red", size=1)+
##    geom_segment(aes(xend=x, yend=0), color="blue", size=0.5)+
##    geom_hline(yintercept=-log10(0.05), linetype="dashed", color="green")+ 
##    xlab("motif")+ylab(bquote(-log[10]~"("~italic(p)~")"))+ 
##    ggtitle(paste(oneMCl, "(lambda.1se)", sep=""))+ 
##    theme_bw()+
##    theme(plot.title=element_text(hjust=0.5))   

## ##
## ###
## figfn <- paste(outdir, "Figure1.4", "_", oneMCl, ".log10p.png", sep="")
## png(figfn, width=480, height=320, res=100)
## print(p2)
## dev.off()


 
###
### bootstrap
###
## fn <- paste(outdir, "1_cell-type_active.peaks.rds", sep="")
## dz <- read_rds(fn)
## npeak <- nrow(dz)

## ###
## ID_bootstrap <- lapply(1:100, function(i){ 
##    set.seed(i)
##    id <- sample(1:npeak, replace=T)
##    id
## })
## ID_bootstrap <- do.call(cbind, ID_bootstrap)

## ### output
## colnames(ID_bootstrap) <- paste("id", 1:100, sep="_")
## opfn2 <- paste(outdir, "3.0_bootstrap.ID.txt", sep="")
## write.table(ID_bootstrap, opfn2, row.names=F, quote=F, sep="\t")



