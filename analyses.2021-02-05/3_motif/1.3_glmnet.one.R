##
rm(list=ls())
###
library(tidyverse)
library(data.table)
library(glmnet)
library(Matrix)
library(doMC)
library(parallel)
library(foreach)
registerDoMC(cores=10)


### parsing arguments
args=commandArgs(trailingOnly=T)
if (length(args)>0){
    ii <- as.numeric(args[1])
    icell <- as.numeric(args[2])
    cell <- as.character(args[3])    
}else{
    ii <- 1
    icell <- 4
    cell <- "Tcell"
}    

###
### a grid of pct values
pct_grid <- c(0.1, 0.05, 0.02, 0.01)
k <- 3
## for (k in 2:length(pct_grid)){
##
pct0 <- pct_grid[k]    
outdir <- paste("./1.3_motif.outs/pct_", pct0, "/", sep="")
if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)

cat(cell, "cv", ii, "pct", pct0, "\n")


### read bootstrap
fn <- "./1.3_motif.outs/3.0_bootstrap.ID.txt"
IDs <- fread(fn, header=T, data.table=F)
IDs <- as.matrix(IDs)
id <- IDs[,ii]


### Y input
fn <- paste("./1.3_motif.outs/pct_", pct0, "/", "1_cell-type_active.peaks.rds", sep="")
dz <- read_rds(fn)

### X input
fn <- "./1.3_motif.outs/0_motif.rds"
motif <- read_rds(fn)


## reshuffle input data
Y <- dz[,icell+1]
Yreshf <- Y[id]
X <- motif[dz$peakAll,]
Xreshf <- X[id,]

## fit lasso
system.time(cvfit <- cv.glmnet(Xreshf, Yreshf, family="binomial", type.measure="class"))
#
b <- coef(cvfit)
b <- as.matrix(b)


### output
opfn <- paste(outdir, "4.", icell, "_", cell, "_Bootstrap_", ii, "_coef_motif.rds", sep="")
write_rds(b, opfn)
##}
### End
