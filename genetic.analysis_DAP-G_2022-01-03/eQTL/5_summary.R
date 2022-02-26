#
library(tidyverse)
library(data.table)
library(clusterProfiler)
library(org.Hs.eg.db)
library(annotables)
library(ggplot2)
library(cowplot)
library(RColorBrewer)

###
###

###
### combine results

outdir <- "./5_summary.outs/dap-g_tss/"
if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)
###
datasets <- read.table("datasets.txt")$V1
geneList <- read.table("zzz_geneList.txt")$V1
###
for (i in 1:length(datasets)){ 
   ##
   condition <- datasets[i]
   cat(i, "condition", condition, "\n") 
   ## 
   resDAP2 <- lapply(geneList, function(ii){
     ###
       fn <- paste("./dap-g_outs/dap-g_tss/", condition, "/", ii, ".SNP.out", sep="")
       if ( file.exists(fn)&file.size(fn)>0){
          res <- fread(fn, header=F, data.table=F)
          res$V8 <- ii
          res0 <- res[1,]
       }else{
          res0 <- NA
       }
       res0
   })
   ###
   resDAP2 <- resDAP2[!is.na(resDAP2)]
   resDAP2 <- do.call(rbind, resDAP2) 
   resDAP2$condition <- condition
   ###
   opfn <- gzfile(paste(outdir, i, "_", condition, "_dap-g.txt.gz", sep=""))
   write.table(resDAP2, opfn, row.names=F, col.names=T, quote=F)
} ###

    

###
### top SNP for each gene
outdir <- "./5_summary.outs/"
datasets <- read.table("datasets.txt")$V1
geneList <- read.table("zzz_geneList.txt")$V1
###
resDAP <- map_dfr(datasets, function(condition){
   ##
   cat(condition, "\n")
   ## 
   resDAP2 <- lapply(geneList, function(ii){
     ###
       fn <- paste("./dap-g_outs/dap-g_tss/", condition, "/", ii, ".SNP.out", sep="")
       if ( file.exists(fn)&file.size(fn)>0){
          res <- fread(fn, header=F, data.table=F)
          res$V8 <- ii
          res0 <- res[1,]
       }else{
          res0 <- NA
       }
       res0
   })
   ###
   resDAP2 <- resDAP2[!is.na(resDAP2)]
   resDAP2 <- do.call(rbind, resDAP2) 
   resDAP2$condition <- condition
   resDAP2
})

res <- resDAP%>%dplyr::arrange(desc(V3))

opfn <- gzfile(paste(outdir, "1_topSNP_dap-g.txt.gz", sep=""))
write.table(resDAP, opfn, row.names=F, col.names=F, quote=F)



###
### reQTL

prefix <-"/wsu/home/groups/piquelab/SCAIP/SCAIP-genetic/" 
###
fn <- paste(prefix, "TWAS_overlap/PTWAS-immune-traits.txt", sep="")
immune <- read.table(fn)$V1
##
fn <- paste(prefix, "TWAS_overlap/media-4_bioRxiv_TableS2.txt", sep="")
twas <- read.table(fn, header=T)%>%mutate(ENSG2=gsub("\\..*", "", Gene))%>%
    filter(Trait%in%immune)


fn <- paste(prefix, "mashr_eQTL/2_reQTL.output/reQTLs_All.infor.rds", sep="")
reQTL <- read_rds(fn)%>%filter(is_TWAS.immune==1)
reQTL2 <- reQTL%>%filter(cell=="Tcell")

###
geneDf <-reQTL2%>%dplyr::select(ENSG2, SYMBOL)%>%distinct(ENSG2,.keep_all=T)



#############################################
### response eGene overlapping twas plots ###
#############################################

pct0 <- 0.02
outdir <- paste("./5_summary.outs/dap-g_pct_", pct0, "/1_example.reGene.outs/", sep="")
if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)


## res <- fread("./5_summary.outs/1_topSNP_dap-g.txt.gz", header=F, data.table=F)
## res <- res%>%arrange(desc(V3))
## res2 <- res%>%dplyr::distinct(V8, .keep_all=T)

## ### top 100 SNPs 
## ## topSNP <- res[1:2000,] 
## anno <- bitr(unique(res2$V8), fromType="ENSEMBL", toType="SYMBOL", OrgDb=org.Hs.eg.db)
## res2 <- res2%>%left_join(anno, by=c("V8"="ENSEMBL"))



### gene
## topSNP2 <- topSNP%>%group_by(V8)%>%top_n(n=1, wt=V3)%>%ungroup()%>%as.data.frame()
## topSNP2 <- topSNP%>%dplyr::distinct(V8,.keep_all=T)
## i <- 500
## ENSG <- topSNP2[i,8]
## SYMBOL <- topSNP2[i,10]
for (i in 1:nrow(geneDf)){

###    
ens <- geneDf$ENSG2[i]
symbol <- geneDf$SYMBOL[i]

cat(ens, symbol, "\n")
    
### data
fig_ls <- lapply(16:20, function(k){
   ###
   condition <- datasets[k] 
   fn <- paste("./dap-g_outs/dap-g_pct_", pct0, "/", condition, "/", ens, ".SNP.out", sep="")

   ###
   if ( file.exists(fn) ){ 
   res <- fread(fn, header=F, data.table=F) 
   res <- res%>%mutate(snpinfor=gsub(";.*", "", V2), rs=gsub(".*;", "", V2))%>%
      separate(snpinfor, into=c("chr", "pos", "ref", "alt"), sep=":")%>%
      mutate(is_sig=ifelse(V3>0.001, 1, 2)) 
   chr <- res$chr[1]
##
   p <- ggplot(res, aes(x=as.numeric(pos), y=V3))+
      geom_point(aes(colour=factor(is_sig), size=factor(is_sig)))+
      scale_size_manual("", values=c("1"=0.8,"2"=0.3))+
      scale_colour_manual("", values=c("1"="red", "2"="grey60"))+
       xlab(bquote("position"~"(chr"~.(chr)~")"))+
      ylab(bquote("PIP"~"("~italic(.(symbol))~")"))+
      ggtitle(bquote(.(condition)))+
      theme(legend.position="none",
         plot.title=element_text(hjust=0.5, size=10),
         panel.background=element_blank(),
         panel.border=element_rect(fill=NA,color="black"),
         panel.grid.major=element_blank(),
         panel.grid.minor=element_blank())
   }else{
       p <- NA
   }    
   p
})
##
##

### print plot    
fig_ls <- fig_ls[!is.na(fig_ls)]
nfig <- length(fig_ls)
figfn <- paste(outdir, "Figure", i, "_", symbol, "_scatter.png", sep="")
###
if (nfig==5){

   png(figfn, width=500, height=480)
   print(plot_grid(plotlist=fig_ls, ncol=2))
   dev.off()
}else if (nfig==4|nfig==3){
   png(figfn, width=500, height=380)
   print(plot_grid(plotlist=fig_ls, ncol=2))
   dev.off()
}else{
   png(figfn, width=500, height=200)
   print(plot_grid(plotlist=fig_ls, ncol=2))
   dev.off()  
}

} ## End




#################################
### summary number of signals ###
#################################

pct0 <- 0.02
outdir <- paste("./5_summary.outs/dap-g_pct_", pct0, "/2_example.multiSignals.outs/", sep="")
if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)

datasets <- read.table("datasets.txt")$V1
geneList <- read.table("zzz_geneList.txt")$V1
###
resDAP <- map_dfr(datasets[16:20], function(condition){
   ##
   cat(condition, "\n")
   ## 
   resDAP2 <- lapply(geneList, function(ii){
     ###
       fn <- paste("./dap-g_outs/dap-g_pct_0.02/", condition, "/", ii, ".cluster.out", sep="")
       if ( file.exists(fn)&file.size(fn)>0){
          res <- read.table(fn, header=F)[,1:4]
          res$ens <- ii
          res$ncluster <- nrow(res)
       }else{
          res <- NA
       }
       res
   })
   ###
   resDAP2 <- resDAP2[!is.na(resDAP2)]
   resDAP2 <- do.call(rbind, resDAP2) 
   resDAP2$condition <- condition
   resDAP2
})


####
#### >2 signals

x <- resDAP%>%filter(ncluster>=4)
gene <- unique(x$ens)
geneDf <- bitr(gene, fromType="ENSEMBL", toType="SYMBOL", OrgDb=org.Hs.eg.db)


for (i in 1:nrow(geneDf)){

###    
ens <- geneDf$ENSEMBL[i]
symbol <- geneDf$SYMBOL[i]

cat(ens, symbol, "\n")
    
### data
fig_ls <- lapply(16:20, function(k){
   ###
   condition <- datasets[k] 
   fn <- paste("./dap-g_outs/dap-g_pct_0.02/", condition, "/", ens, ".SNP.out", sep="")

   ###
   if ( file.exists(fn) ){ 
   res <- fread(fn, header=F, data.table=F) 
   res <- res%>%mutate(snpinfor=gsub(";.*", "", V2), rs=gsub(".*;", "", V2))%>%
      separate(snpinfor, into=c("chr", "pos", "ref", "alt"), sep=":")%>%
      mutate(is_sig=ifelse(V3>0.001, 1, 2)) 
   chr <- res$chr[1]
##
   p <- ggplot(res, aes(x=as.numeric(pos), y=V3))+
      geom_point(aes(colour=factor(is_sig), size=factor(is_sig)))+
      scale_size_manual("", values=c("1"=0.8,"2"=0.3))+
      scale_colour_manual("", values=c("1"="red", "2"="grey60"))+
       xlab(bquote("position"~"(chr"~.(chr)~")"))+
      ylab(bquote("PIP"~"("~italic(.(symbol))~")"))+
      ggtitle(bquote(.(condition)))+
      theme(legend.position="none",
         plot.title=element_text(hjust=0.5, size=10),
         panel.background=element_blank(),
         panel.border=element_rect(fill=NA,color="black"),
         panel.grid.major=element_blank(),
         panel.grid.minor=element_blank())
   }else{
       p <- NA
   }    
   p
})
##
##

### print plot    
fig_ls <- fig_ls[!is.na(fig_ls)]
nfig <- length(fig_ls)
figfn <- paste(outdir, "Figure", i, "_", symbol, "_scatter.png", sep="")
###
if (nfig==5){

   png(figfn, width=500, height=480)
   print(plot_grid(plotlist=fig_ls, ncol=2))
   dev.off()
}else if (nfig==4|nfig==3){
   png(figfn, width=500, height=380)
   print(plot_grid(plotlist=fig_ls, ncol=2))
   dev.off()
}else{
   png(figfn, width=500, height=200)
   print(plot_grid(plotlist=fig_ls, ncol=2))
   dev.off()  
}

} ## End



##################
### locus zoom ###         
##################

outdir <- "./5_summary.outs/dap-g_pct_0.02/3_example.zoom.outs/"
if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)


x <- resDAP%>%filter(ncluster>=2)
gene <- unique(x$ens)
geneDf <- bitr(gene, fromType="ENSEMBL", toType="SYMBOL", OrgDb=org.Hs.eg.db)

## "SLC16A1"
## "IL7R"
## "HLA-B"
## "HLA-C"
## "S100B"
symbol <- "SLC16A1"
x <- geneDf%>%filter(SYMBOL==symbol)
ens <- x[1,1]


### read data
res <- map_dfr(16:20, function(i){
###    
condition <- datasets[i] 
fn <- paste("./dap-g_outs/dap-g_pct_0.02/", condition, "/", ens, ".SNP.out", sep="")
if ( file.exists(fn)&file.size(fn)>0 ){ 
   res <- fread(fn, header=F, data.table=F) 
   res <- res%>%mutate(snpinfor=gsub(";.*", "", V2), rs=gsub(".*;", "", V2))%>%
      separate(snpinfor, into=c("chr", "pos", "ref", "alt"), sep=":")%>%
      mutate(is_sig=ifelse(V3>-1, 1, 0),"conditions"=condition) 
}
res <- res%>%mutate(pos=as.numeric(pos))
res
})
##
x <- res%>%filter(V5!=-1)
xmin <- min(x$pos)-(0.1e+05)
xmax <- max(x$pos)+(0.1e+05)
res <- res%>%filter(pos>=xmin, pos<=xmax)
res <- res%>%mutate(pos2=pos/(1e+06))

###
### plots for 5 conditions
fig_ls <- lapply(16:20, function(i){
###    
condition <- datasets[i] 

res2 <- res%>%filter(conditions==condition)
chr <- res2$chr[1]    
##
p <- ggplot(res2, aes(x=as.numeric(pos2), y=V3))+
   geom_point(aes(colour=factor(V5), size=factor(is_sig)))+
   scale_size_manual("", values=c("1"=1,"0"=0.3))+
   scale_colour_manual("", values=c("-1"="grey60","1"="red", "2"="green",
       "3"="blue", "4"="#984ea3", "5"="#ffaa00"))+
   guides(size="none")+ 
   xlab(bquote("position"~"(chr"~.(chr)~")"))+
   ylab(bquote("PIP"~"("~italic(.(symbol))~")"))+
   ggtitle(bquote(.(condition)))+
   theme(
       ## legend.background=element_blank(),
       ## legend.box.background=element_blank(),
       ## legend.key=element_blank(),
       legend.position="none",
       plot.title=element_text(hjust=0.5, size=12),
       panel.background=element_blank(),
       panel.border=element_rect(fill=NA,color="black"),
       panel.grid.major=element_blank(),
       panel.grid.minor=element_blank())
p
})


###
###one example
p <- ggplot(res%>%filter(conditions=="Tcell_CTRL"), aes(x=pos2, y=V3))+
   geom_point(aes(colour=factor(V5), size=factor(is_sig)))+
   scale_size_manual("", values=c("1"=1,"0"=0.3))+
   scale_colour_manual("", values=c("-1"="grey60","1"="red", "2"="green",
       "3"="blue", "4"="#984ea3", "5"="#ffaa00"))+
   guides(size="none", colour=guide_legend(ncol=2, override.aes=list(size=2.5)))+ 
   xlab(bquote("position"~"(chr"~.(chr)~")"))+
   ylab(bquote("PIP"~"("~italic(.(symbol))~")"))+
   ggtitle(bquote(.(condition)))+
   theme(legend.key=element_blank(),
       plot.title=element_text(hjust=0.5, size=12),
       panel.background=element_blank(),
       panel.border=element_rect(fill=NA,color="black"),
       panel.grid.major=element_blank(),
       panel.grid.minor=element_blank())

fig_ls[[6]] <- get_legend(p)
    
figfn <- paste(outdir, "Figure_", symbol, "_scatter.png", sep="")       
png(figfn, width=500, height=480)
print(plot_grid(plotlist=fig_ls, ncol=2))
dev.off()
### End


