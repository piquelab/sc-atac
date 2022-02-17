##
library(tidyverse)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(viridis)


####################################################################
### forste plot, shared the same annotation across 20 conditions ###
####################################################################

outdir <- "./3_summary.outs/torus_peak_condition/"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)


datasets <- read.table("datasets.txt")
est2 <- map_dfr(1:nrow(datasets), function(i){
###
  condition <- datasets$V1[i]
   fn <- paste("./torus_output/torus_peak_condition/", condition, ".est", sep="")
   est <- read.table(fn)
   est <- est[2,2:4]
   est <- est%>%mutate(conditions=condition,
      MCls=gsub("_.*", "", conditions), treats=gsub(".*_|-EtOH", "", condition))    
###
})
est2 <- est2%>%mutate(conditions2=gsub("-", "+", gsub("-EtOH", "", conditions))) 
####   
###    
p <- ggplot(est2, aes(x=V2, y=conditions2, color=factor(MCls)))+
   geom_errorbarh(aes(xmax=V4, xmin=V3), size=0.5, height=0.2)+
   geom_point(shape=19, size=0.5)+
   geom_vline(aes(xintercept=0), size=0.25, linetype="dashed")+
   xlab("log odds ratio")+
   ## scale_y_discrete(labels=ylabel.names)+
   scale_colour_manual(
      values=c("Bcell"="#4daf4a", "Monocyte"="#984ea3",
               "NKcell"="#aa4b56", "Tcell"="#ffaa00"))+
   ggtitle("Enrichment of peaks")+  
   theme_bw()+
   theme(plot.title=element_text(hjust=0.5, size=10),
         legend.position="none",
         axis.title.x=element_text(size=10),
         axis.title.y=element_blank(),
         axis.text=element_text(size=8))
###
figfn <- paste(outdir, "Figure1.1_est.png", sep="")
png(figfn, width=360, height=500, res=120)
print(p)
dev.off()





########################################################################
### forest plot for sharing the same annotation across 20 conditions ###
########################################################################


### 0 representing SNPs are not at any motif
### 1 motif
### 2 representing SNPs in the motif and peaks

outdir <- "./3_summary.outs/torus_motif3_all2/"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)


datasets <- read.table("datasets.txt")
est2 <- map_dfr(1:nrow(datasets), function(i){
###
   condition <- datasets$V1[i]
   fn <- paste("./torus_output/torus_motif3_all2/", condition, ".est", sep="")
   est <- read.table(fn)
   est <- est[2:3,1:4]
   est <- est%>%mutate(conditions=condition,
      MCls=gsub("_.*", "", conditions), treats=gsub(".*_|-EtOH", "", condition))    
###
})
est2 <- est2%>%mutate(conditions2=gsub("-", "+", gsub("-EtOH", "", conditions))) 
####   
###    
p <- ggplot(est2, aes(x=V2, y=conditions2, color=factor(MCls), alpha=factor(V1)))+
   geom_errorbarh(aes(xmax=V4, xmin=V3), size=0.5, height=0.2)+
   geom_point(shape=19, size=0.5)+
   geom_vline(aes(xintercept=0), size=0.25, linetype="dashed")+
   xlab("log odds ratio")+
   ## scale_y_discrete(labels=ylabel.names)+
   scale_colour_manual(
      values=c("Bcell"="#4daf4a", "Monocyte"="#984ea3",
               "NKcell"="#aa4b56", "Tcell"="#ffaa00"))+
   scale_alpha_manual(values=c("binding.1"=0.4, "binding.2"=1),
       labels=c("binding.1"="motif", "binding.2"="motif&peak"))+
   guides(colour="none", alpha=guide_legend(title=NULL, override.aes=list(colour="#984ea3")))+ 
   ggtitle("Enrichment of motifs")+  
   theme_bw()+
   theme(plot.title=element_text(hjust=0.5, size=10),
         axis.title.x=element_text(size=10),
         axis.title.y=element_blank(),
         axis.text=element_text(size=8))
###
figfn <- paste(outdir, "Figure1.1_est.png", sep="")
png(figfn, width=550, height=500, res=120)
print(p)
dev.off()



### 0 representing SNPs are not at any peaks
### 1 peaks but not motifs
### 2 representing SNPs in the peaks and motif

outdir <- "./3_summary.outs/torus_peak_all2/"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)


datasets <- read.table("datasets.txt")
est2 <- map_dfr(1:nrow(datasets), function(i){
###
   condition <- datasets$V1[i]
   fn <- paste("./torus_output/torus_peak_all2/", condition, ".est", sep="")
   est <- read.table(fn)
   est <- est[2:3,1:4]
   est <- est%>%mutate(conditions=condition,
      MCls=gsub("_.*", "", conditions), treats=gsub(".*_|-EtOH", "", condition))    
###
})
est2 <- est2%>%mutate(conditions2=gsub("-", "+", gsub("-EtOH", "", conditions))) 
####   
###
est3 <- est2%>%mutate(V3=ifelse(V3<(-6), -6, V3), V4=ifelse(V4>6, 6, V4))
p <- ggplot(est3, aes(x=V2, y=conditions2, color=factor(MCls), alpha=factor(V1)))+
   geom_errorbarh(aes(xmax=V4, xmin=V3), size=0.5, height=0.2)+
   geom_point(shape=19, size=0.5)+
   geom_vline(aes(xintercept=0), size=0.25, linetype="dashed")+
   xlab("log odds ratio")+xlim(-6,6)+
   ## scale_y_discrete(labels=ylabel.names)+
   scale_colour_manual(
      values=c("Bcell"="#4daf4a", "Monocyte"="#984ea3",
               "NKcell"="#aa4b56", "Tcell"="#ffaa00"))+
   scale_alpha_manual(values=c("peaking.1"=0.4, "peaking.2"=1),
       labels=c("peaking.1"="peak", "peaking.2"="peak&motif"))+
   guides(colour="none", alpha=guide_legend(title=NULL, override.aes=list(colour="#984ea3")))+ 
   ggtitle("Enrichment of peak")+  
   theme_bw()+
   theme(plot.title=element_text(hjust=0.5, size=10),
         axis.title.x=element_text(size=10),
         axis.title.y=element_blank(),
         axis.text=element_text(size=8))
###
figfn <- paste(outdir, "Figure1.1_est.png", sep="")
png(figfn, width=550, height=500, res=120)
print(p)
dev.off()




###########################
#### Peaks by cell type ###
###########################



outdir <- "./3_summary.outs/torus_peak_cell/"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)


datasets <- read.table("datasets.txt")

##
est <- map_dfr (c("Bcell", "Monocyte", "NKcell", "Tcell"), function(cell){
   ###  
   cat(cell, "\n")
   est2 <- map_dfr(1:nrow(datasets), function(i){
      condition <- datasets$V1[i]
      fn <- paste("./torus_output/torus_peak_cell/", cell, "/", condition, ".est", sep="")
      est0 <- read.table(fn)
      est0 <- est0[2,1:4]
      est0 <- est0%>%mutate(conditions=condition,
          MCls=gsub("_.*", "", conditions), treats=gsub(".*_|-EtOH", "", condition), peak2=cell)    
###
      est0
   })
   est2 <- est2%>%mutate(conditions2=gsub("-", "+", gsub("-EtOH", "", conditions))) 
})


####   
###    
p <- ggplot(est, aes(x=V2, y=conditions2, color=factor(MCls)))+
   geom_errorbarh(aes(xmax=V4, xmin=V3), size=0.5, height=0.2)+
   geom_point(shape=19, size=0.5)+
   geom_vline(aes(xintercept=0), size=0.25, linetype="dashed")+
   xlab("log odds ratio")+xlim(-5,5)+
   ## scale_y_discrete(labels=ylabel.names)+
   scale_colour_manual(
      values=c("Bcell"="#4daf4a", "Monocyte"="#984ea3",
               "NKcell"="#aa4b56", "Tcell"="#ffaa00"))+
   facet_wrap(~peak2, ncol=4)+
   ggtitle("Enrichment of peaks")+  
   theme_bw()+
   theme(plot.title=element_text(hjust=0.5, size=10),
         legend.position="none",
         axis.title.x=element_text(size=10),
         axis.title.y=element_blank(),
         axis.text=element_text(size=8))
###
figfn <- paste(outdir, "Figure1.1_est.png", sep="")
png(figfn, width=850, height=420, res=120)
print(p)
dev.off()



#################################################################
### forest plot for cell-type and response motids annotation  ###
#################################################################

#### final version

rm(list=ls())

pct_grid <- c(0.05, 0.02, 0.01)[2]

for (pct0 in pct_grid){
    
option <- paste("torus_peak_pct_", pct0, "_response3", sep="")
## option <- "torus_peak_condition"
    
###
outdir <- paste("./3_summary.outs/", option, "/", sep="")
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)

####
datasets <- read.table("datasets.txt")
est2 <- map_dfr(16:20, function(i){
###
  condition <- datasets$V1[i]
  fn <- paste("./torus_output/", option, "/", condition, ".est", sep="") 
  est <- read.table(fn)
  est <- est%>%mutate(order=1:nrow(est), conditions=condition,
      MCls=gsub("_.*", "", conditions),
      treats=gsub(".*_|-EtOH", "", conditions),
      label.name=gsub("]", "", gsub(".*\\[", "", V1)))
   est <- est[c(2:4,8:21),]  
})


###
tmp <- est2%>%filter(conditions=="Tcell_CTRL")
## ylabel.names <- as.character(tmp$label.name)
## names(ylabel.names) <- tmp$order

ylab.name <- gsub("]", "", gsub(".*\\[", "", tmp$V1))
ylab.name[1:3] <- c("peak_nomotif","peak&cell-type-motif", "peak&condition-motif") 
names(ylab.name) <- as.character(tmp$order)


###
###
ii <- "Tcell"
p <- ggplot(est2%>%filter(MCls==ii), aes(x=V2, y=as.factor(order),color=factor(treats)))+
      geom_errorbarh(aes(xmax=V4, xmin=V3), size=0.5, height=0.2)+
      geom_point(shape=19, size=0.5)+
      geom_vline(aes(xintercept=0), size=0.25, linetype="dashed")+
      scale_x_continuous("log odds ratio", breaks=seq(-8, 6, 2), limits=c(-8.5,6.5)) +
      scale_y_discrete("SNP annotation", labels=ylab.name)+
      scale_colour_manual(
          values=c("CTRL"="#828282", "LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
                  "PHA"="#a6cee3", "PHA-DEX"="#1f78b4"))+
      ggtitle(ii)+
      facet_wrap(~treats, ncol=5)+
      theme_bw()+
      theme(plot.title=element_text(hjust=0.5, size=10),
            legend.position="none",
            axis.title=element_text(size=10),
            axis.text=element_text(size=8))
###
figfn <- paste(outdir, "Figure0.1_", ii, ".tss.est.png", sep="")
png(figfn, width=900, height=500, res=120)
print(p)
dev.off()

} ## End    
