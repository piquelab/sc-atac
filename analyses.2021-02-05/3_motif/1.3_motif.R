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


outdir <- "./1.3_motif.outs/"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)


###############################
### cell-type active peaks  ###
###############################

## rm(list=ls())

### motif
fn <- "./1_motif.outs/1_scATAC.motif.rds" 
atac <- read_rds(fn)
motif <- Motifs(atac)
X <- motif@data
##
opfn <- paste(outdir, "0_motif.rds", sep="")
write_rds(X, opfn)

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
MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
###
dz <- map_dfc(1:length(MCls), function(i){
###
   oneMCl <- MCls[i]
   atac2 <- subset(atac, subset=MCls==oneMCl)
   count2 <- atac2@assays$ATAC@counts
   rpz <- rowMeans(count2>0)
###
   peakSel <- rownames(count2)[rpz>0.1]
   z <- ifelse(peakAll%in%peakSel, 1, 0)
   z
})
dz <- cbind(peakAll, dz)
names(dz)[2:5] <- paste("d_", MCls, sep="")
rownames(dz) <- dz$peakAll
###
opfn <- paste(outdir, "1_cell-type_active.peaks.rds", sep="")
write_rds(dz, opfn)



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
outdir <- "./1.3_motif.outs/"

fn <- paste(outdir, "1_cell-type_active.peaks.rds", sep="")
dz <- read_rds(fn)

###
fn <- paste(outdir, "0_motif.rds", sep="")
motif <- read_rds(fn)


MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")

###
### loop by cell-type
B <- lapply(1:4, function(i){
###
   oneMCl <- MCls[i]
   cat(oneMCl, i, "\n")
###    
   Y <- dz[,i+1]
   X <- motif[dz$peakAll,]
## fit lasso
   system.time(cvfit <- cv.glmnet(X, Y, family="binomial", type.measure="class"))
#
   opfn <- paste(outdir, "2_", i, "_", oneMCl, ".glmnet.rds", sep="")
   write_rds(cvfit, opfn)

##    
   b <- coef(cvfit, s="lambda.min")
   b <- as.matrix(b)
})
B <- do.call(cbind, B)

###
colnames(B) <- paste("coef_", MCls, sep="")
##
opfn <- paste(outdir, "2_coef_motif.rds", sep="")
write_rds(B, opfn)


##
## cvplot
fn <- paste(outdir, "2_4_Tcell.glmnet.rds", sep="")
cvfit <- read_rds(fn)
###
figfn <- paste(outdir, "Figure0.4_Tcell.lambda.png", sep="")
png(figfn, width=520, height=380, res=100)
plot(cvfit)
dev.off()


###
### coef plot
fn <- paste(outdir, "2_coef_motif.rds", sep="")
B <- read_rds(fn)
B2 <- B[-1,]

MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
i <- 4
oneMCl <- MCls[i]
##

df <- data.frame(x=1:nrow(B2), y=B2[,i])
p <- ggplot(df, aes(x,y))+
   geom_point(color="red", size=1)+
   geom_segment(aes(xend=x, yend=0), color="blue", size=0.5)+
   xlab("motif")+ylab("coefficients")+ 
   ggtitle(oneMCl)+ 
   theme_bw()+
   theme(plot.title=element_text(hjust=0.5))

###
figfn <- paste(outdir, "Figure1.", i, "_", oneMCl, ".coef.png", sep="")
png(figfn, width=480, height=420, res=120)
print(p)
dev.off()


###
motif.name <- unlist(motif@motif.names)
dfmotif <- data.frame(motif.id=as.character(rownames(B2)), coef=B2[,4])%>%
    mutate("motif.name"=motif.name[motif.id]) 
x <- dfmotif%>%arrange(desc(coef))
                       
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



