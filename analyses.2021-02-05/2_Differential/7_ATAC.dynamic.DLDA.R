###
library(Matrix)
library(MASS)
library(scales)
library(tidyverse)
library(parallel)
library(data.table)
##
library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(Signac)
#library(harmony)
library(annotables) 
library(org.Hs.eg.db)
###
library(ggplot2)
library(cowplot, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
library(grid)
library(gridExtra)
library(ggExtra)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(viridis)
library(ggrastr)
theme_set(theme_grey())

rm(list=ls())

outdir <- "./7_dynamic.outs/"
if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=FALSE)


############################################################
###      psedotime definition based on 6,571 DEGs        ###
###       diagonal linear Discriminant analysis          ###
###       Last modified by Julong wei, 2022-04-15        ###
############################################################



########################################################
###  extract atac data for each cell-type separately ### 
########################################################

### differential peaks
res <- read_rds("./1.2_DiffPeak.outs/2.0_DESeq.results.rds")%>%as.data.frame()
res <- res%>%mutate(gamma=ifelse((abs(estimate)>0.5)&(p.adjusted<0.1), 1, 0))

peakSel <- res%>%filter(gamma==1)%>%dplyr::pull(gene)%>%unique()


### read seurat object
atac <- read_rds("../1_processing/5.1_reCallPeak.outs/3_scATAC.annot.rds")
atac2 <- subset(atac, features=peakSel)
 
atac2 <- SplitObject(atac2, split.by="MCls")

###
### normalize data for each cell-type separately
MCls <- sort(unique(atac$MCls))[-2]
for (oneMCl in MCls){
   ##
   cat(oneMCl, "\n") 
   atac0 <- atac2[[oneMCl]]
   atac0 <- atac0%>%
       RunTFIDF()%>%
       FindTopFeatures(min.cutoff="q0")%>%
       RunSVD(n=100)%>%
       RunUMAP(reduction="lsi", dims=2:50, reduction.name="umap.atac", reduction.key="atacUMAP")
   ###
   opfn <- paste(outdir, "1_MCl.", oneMCl, ".cluster.rds", sep="")
   write_rds(atac0, opfn) 
}



######################
### Calculate DLDA ###
######################

### Differential peaks

res <- read_rds("./1.2_DiffPeak.outs/2.0_DESeq.results.rds")%>%as.data.frame()
res <- res%>%dplyr::rename(beta=estimate)%>%
    mutate(gamma=ifelse((abs(beta)>0.5)&(p.adjusted<0.1), 1, 0))

MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
option <- "DiagLDA2"

outdir2 <- paste("./7_dynamic.outs/", option, sep="")
if ( !file.exists(outdir2)) dir.create(outdir2, showWarnings=F, recursive=FALSE)

tmp <- lapply(1:3, function(i){
###
   oneMCl <- MCls[i]
   s1 <- Sys.time()
   fn <- paste("./7_dynamic.outs/1_MCl.", oneMCl, ".cluster.rds", sep="")
   sc <- read_rds(fn)
   X <- sc@assays$ATAC@data
   rn <- rownames(X)
    
   meta <- sc@meta.data%>%
       mutate(treats=gsub(".*-ATAC-|_.*", "", NEW_BARCODE))  
   treats <- unique(meta$treats)
   ntreat <- length(treats)
   ncell <- ncol(X)
   np <- nrow(X) 
    
   ### variance across genes  
   ## Xc <- map_dfc(treats, function(one){
   ##    ii <- as.character(meta[meta$treats==one, "NEW_BARCODE"])
   ##    Xi <- X[,ii]
   ##    mu_k <- apply(Xi, 1, mean)
   ##    xc <- sweep(Xi, 1, mu_k, "-")     
   ##    rowSums(xc*xc)      
   ## })
   mu <- rowMeans(X)
   X0 <- matrix(mu, np, 1)%*%matrix(1, 1, ncell)

   ##
   index <- list()
   step <- trunc(ncell/10) 
   for ( ii in 1:10){
      i1 <- (ii-1)*step+1
      i2 <- ifelse(ii==10, ncell, (ii-1)*step+step)
      index[[ii]] <- c(i1, i2)
   } ##
    
   splitMxList <- lapply(1:10, function(ii){
       ##
       i1 <- index[[ii]][1]
       i2 <- index[[ii]][2]
       mx <- X[,i1:i2]-X0[,i1:i2]
       mx
    })

   Xc <- do.call(cbind, splitMxList)
    
   tmp <- Xc*Xc
    
   var <- rowSums(tmp)/(ncell-1)
   names(var) <- rn 

   ###
   contrast_ls <- list("LPS"=c("CTRL", "LPS"),
                       "LPS-DEX"=c("LPS", "LPS-DEX"),
                       "PHA"=c("CTRL", "PHA"),
                       "PHA-DEX"=c("PHA", "PHA-DEX"))
   contrast_nn <- names(contrast_ls)
   metaNew <- sapply(contrast_nn, function(nn){
      ### 
      oneX <- contrast_ls[[nn]] 
      ### extract lfc as weights
      res0 <- res%>%filter(MCls==oneMCl, contrast==nn)
      ## gamma <- rep(0, length(rn))
      ## names(gamma) <- rn
      ## gamma[res0$gene] <- 1

      ## res0 <- res%>%filter(MCls==oneMCl, contrast==nn)
      ## lfc <- abs(res0$beta)
      ## lfc[is.na(lfc)] <- 0
      ## names(lfc) <- as.character(res0$gene)
      ## lfc0 <- lfc[rn]
      
      ### matrix 1 
      ii <- meta%>%filter(treats==oneX[1])%>%dplyr::pull(NEW_BARCODE)
      X1 <-X[,ii]
      mu1 <- rowMeans(X1)

      ### matrix 2
      ii <- meta%>%filter(treats==oneX[2])%>%dplyr::pull(NEW_BARCODE)
      X2 <-X[,ii]
      mu2 <- rowMeans(X2)


      ###
      if ( option=="DiagLDA"){
         s0 <- rep(1, length(rn))
      }
      
      ###
      if( option=="DiagLDA2"){
         ## 
         if ( sum(res0$gamma>0, na.rm=T)>0){ 
            gamma <- res0$gamma
            names(gamma) <- res0$gene
            gamma[is.na(gamma)] <- 0
            s0 <- gamma[rn]
         }else{
           s0 <- rep(1, length(rn))  
         }
         ### 
      }
      
      ###
      if (option=="DiagLDA3"){
         lfc <- abs(res0$beta)
         names(lfc) <- res0$gene
         lfc[is.na(lfc)] <- 0
         s0 <- lfc[rn]
      }  
          

      Diff <- as.matrix((mu2-mu1)*(1/var)*s0)
      Diff[is.na(Diff)] <- 0
      ##
      ## X12 <- cbind(X1, X2)
      ## mu12 <- apply(X12, 1, mean)
      ## X12 <- sweep(X12, 1, mu12, "-")
      z1 <- as.vector(crossprod(Xc, Diff))
      z1
   })
    
   metaNew <- as.data.frame(metaNew) 
   names(metaNew) <- paste("z_", contrast_nn, sep="")  
   metaNew <- metaNew%>%
      mutate(NEW_BARCODE=colnames(X),
             treats=gsub(".*-ATAC-|_.*", "", NEW_BARCODE),
             treat2=gsub("-EtOH", "", treats)) 
 
   ###
   opfn <- paste(outdir, option, "/2_MCls.", oneMCl, ".DLDA.rds", sep="")
   write_rds(metaNew, opfn)
      
   s2 <- Sys.time()
   d12 <- difftime(s2, s1, units="mins")
   cat(oneMCl, ":", d12, "\n")
})



#########################
### summary LDA plots ###
#########################
### scatter plots plus density plots


