##
library(tidyverse)
library(data.table)
library(qvalue)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(RColorBrewer)
library(viridis)


### For example Tcell_LPS
### Tcell peaks, Tcell motifs, differential motifs in LPS across cell-types

 
outdir <- "./3_summary.outs/torus_peak_combineNew/"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)


###########################################################################
### 1. peaks, cell-type motifs, each treatment motifs across cell-types ###
###########################################################################

fn <- "./torus_output/torus_peak_combineNew/ALOFT_combine.est"
est <- read.table(fn)
est2 <- est[2:10,]
est2 <- est2%>%mutate(order=as.numeric(1:nrow(est2)), category=as.character(gsub(".1", "", V1)))%>%
   dplyr::rename("odds"="V2", "CI_lower"="V3", "CI_upper"="V4")
est2 <- est2%>%mutate(category2=fct_reorder(category, order))


col2 <- c("Bcell"="#4daf4a", "Monocyte"="#984ea3","NKcell"="#aa4b56", "Tcell"="#ffaa00",
  "LPS"="#fb9a99", "LPS-DEX"="#e31a1c", "PHA"="#a6cee3", "PHA-DEX"="#1f78b4", "peaking"="#e7298a")

###
p <- ggplot(est2, aes(x=odds, y=category2, color=factor(category2)))+
    geom_errorbarh(aes(xmax=CI_upper, xmin=CI_lower), size=0.5, height=0.2)+
    geom_point(shape=19, size=0.5)+
    geom_vline(aes(xintercept=0), size=0.25, linetype="dashed")+
##    scale_y_discrete(labels=ylab)+
    scale_x_continuous("log odds ratio", limits=c(-1, 3))+
    scale_colour_manual(values=col2)+
    theme_bw()+
    theme(legend.title=element_blank(),
          legend.key.size=grid::unit(0.8, "lines"),
          axis.title.y=element_blank(),
          axis.title.x=element_text(size=10),
          strip.text=element_text(size=12),
          ##plot.title=element_text(hjust=0.5, size=12),
          ## axis.text.x=element_text(size=9, angle=-90, hjust=0, vjust=0.5),
          axis.text=element_text(size=9))
 
figfn <- paste(outdir, "Figure1.1_combine_annot.est.png", sep="")
png(figfn, width=480, height=380, res=120)
p
dev.off()




################################################################################################
### 2. combine annotation, peaks, cell-type motifs, treatment motifs in different cell-types ###
################################################################################################

fn <- "./torus_output/torus_peak_combineNew/ALOFT_combine2.est"
est <- read.table(fn)
est2 <- est[2:22,]
est2 <- est2%>%mutate(order=as.numeric(1:nrow(est2)),
   category=as.character(gsub(".1", "", V1)),
   gr=ifelse(grepl("_", category), gsub(".*_", "", category), category))%>%
   dplyr::rename("odds"="V2", "CI_lower"="V3", "CI_upper"="V4")
est2 <- est2%>%mutate(category2=fct_reorder(category, order))


col2 <- c("Bcell"="#4daf4a", "Monocyte"="#984ea3","NKcell"="#aa4b56", "Tcell"="#ffaa00",
  "LPS"="#fb9a99", "LPS-DEX"="#e31a1c", "PHA"="#a6cee3", "PHA-DEX"="#1f78b4", "peaking"="#e7298a")

p2 <- ggplot(est2, aes(x=odds, y=category2, color=factor(gr)))+
    geom_errorbarh(aes(xmax=CI_upper, xmin=CI_lower), size=0.5, height=0.2)+
    geom_point(shape=19, size=0.5)+
    geom_vline(aes(xintercept=0), size=0.25, linetype="dashed")+
    ##scale_y_discrete(labels=ylab)+
    scale_x_continuous("log odds ratio", limits=c(-1, 3))+
    scale_colour_manual(values=col2)+
    theme_bw()+
    theme(legend.title=element_blank(),
          legend.key.size=grid::unit(0.8, "lines"),
          axis.title.y=element_blank())
          ##axis.title.x=element_text(size=10),
          ###strip.text=element_text(size=12),
          ##plot.title=element_text(hjust=0.5, size=12),
          ## axis.text.x=element_text(size=9, angle=-90, hjust=0, vjust=0.5),
          ###axis.text=element_text(size=9))
 
figfn <- paste(outdir, "Figure1.2_combine2_annot.est.png", sep="")
png(figfn, width=520, height=500, res=120)
p2
dev.off()


######################################################
### 3. peaks, cell-type motifs and response motifs ###
######################################################

fn <- "./torus_output/torus_peak_combineNew/ALOFT_Union.est"
est <- read.table(fn)
est2 <- est[2:4,]
est2 <- est2%>%mutate(order=as.numeric(1:nrow(est2)))%>%
   dplyr::rename("peaking_d"="V1", "odds"="V2", "CI_lower"="V3", "CI_upper"="V4")

p3 <- ggplot(est2, aes(x=odds, y=as.character(order), color=factor(order)))+
    geom_errorbarh(aes(xmax=CI_upper, xmin=CI_lower), size=0.5, height=0.2)+
    geom_point(shape=19, size=0.5)+
    geom_vline(aes(xintercept=0), size=0.25, linetype="dashed")+
    scale_y_discrete(labels=c("1"="Peaks", "2"="Cell-type motifs", "3"="Treatment motifs"))+
    scale_x_continuous("log odds ratio", limits=c(-1, 3))+
    scale_colour_manual(values=c("1"="#7570b3", "2"="#1b9e77", "3"="#e7298a"),
       labels=c("1"="Peaks", "2"="Cell-type motifs", "3"="Treatment motifs"))+
    theme_bw()+
    theme(legend.title=element_blank(),
          legend.key.size=grid::unit(0.8, "lines"),
          axis.title.y=element_blank())
          ##axis.title.x=element_text(size=10),
          ###strip.text=element_text(size=12),
          ##plot.title=element_text(hjust=0.5, size=12),
          ## axis.text.x=element_text(size=9, angle=-90, hjust=0, vjust=0.5),
          ###axis.text=element_text(size=9))
 
figfn <- paste(outdir, "Figure1.3_Union_annot.est.png", sep="")
png(figfn, width=480, height=300, res=120)
p3
dev.off()



##########################
### meta T-cells scaip ###
##########################
fn <- "/nfs/rprdata/julong/sc-atac/genetic.analysis_torus_2021-10-14/eQTL/datasets.txt" 
datasets <- read.table(fn)$V1
est3 <- map_dfr(16:20, function(i){
    ##
    condition <- datasets[i]
    fn <- paste("/nfs/rprdata/julong/sc-atac/genetic.analysis_torus_2021-10-14/eQTL/torus_output/torus_peak_union/", condition, ".est", sep="")
    est <- read.table(fn)
    est <- est[2:4,]
    est <- est%>%mutate(conditions=condition,
       MCls=gsub("_.*", "", conditions),
       order=gsub("peaking.", "", V1))
###
    est
})
###
est3 <- est3%>%dplyr::rename("odds"="V2", "CI_lower"="V3", "CI_upper"="V4")%>%
    mutate(se=abs(odds-CI_lower)/1.96)

DF <- NULL
for (i in 1:3){
   ##
   est0 <- est3%>%dplyr::filter(order==i)    
   b <- est0$odds 
   v <- (est0$se)^2
   ##
   bmeta <- sum((1/v)*(1/sum(1/v))*b)
   vmeta <- 1/sum(1/v)
   se <- sqrt(vmeta) 
   df2 <- data.frame(order=i, b=bmeta, se=sqrt(vmeta), CI_lower=bmeta-1.96*se, CI_upper=bmeta+1.96*se)
   DF <- rbind(DF, df2)
}

DF$MCls <- "Tcell"
DF <- DF%>%dplyr::select(odds=b, CI_lower, CI_upper, order, MCls)

####
#### ALOFT
fn <- "./torus_output/torus_peak_combineNew/ALOFT_Union.est"
est <- read.table(fn)
est2 <- est[2:4,]
est2 <- est2%>%mutate(order=as.numeric(1:nrow(est2)), MCls="Bulk")%>%
   dplyr::rename("peaking_d"="V1", "odds"="V2", "CI_lower"="V3", "CI_upper"="V4")%>%
   dplyr::select(odds, CI_lower, CI_upper, order, MCls)

plotDF <- rbind(DF, est2)%>%
    mutate(comb=paste(order, MCls, sep="_"))

###
###
p <- ggplot(plotDF,
    aes(x=odds, y=as.character(order), color=factor(order), group=factor(MCls), alpha=factor(MCls)))+
    geom_errorbarh(aes(xmax=CI_upper, xmin=CI_lower), size=0.5, height=0.2, position=position_dodge(0.5))+
    geom_point(shape=19, size=0.5, position=position_dodge(0.5))+
    geom_vline(aes(xintercept=0), size=0.25, linetype="dashed")+
    scale_y_discrete(labels=c("1"="Peaks", "2"="Cell-type motifs", "3"="Treatment motifs"))+
    scale_x_continuous("log odds ratio", limits=c(-1, 4.5))+
    scale_colour_manual("", values=c("1"="#7570b3", "2"="#1b9e77", "3"="#e7298a"),
       labels=c("1"="Peaks", "2"="Cell-type motifs", "3"="Treatment motifs"), guide="none")+
   scale_alpha_manual("", values=c("Tcell"=1, "Bulk"=0.5),
       labels=c("Tcell"="T-cells (SCAIP)", "Bulk"="Whole blood (ALOFT)"),               
       guide=guide_legend(override.aes=list(color="#e7298a")))+ 
    theme_bw()+
    theme(legend.title=element_blank(),
          legend.key.size=grid::unit(1, "lines"),
          axis.title.y=element_blank())
          ##axis.title.x=element_text(size=10),
          ###strip.text=element_text(size=12),
          ##plot.title=element_text(hjust=0.5, size=12),
          ## axis.text.x=element_text(size=9, angle=-90, hjust=0, vjust=0.5),
          ###axis.text=element_text(size=9))
 
figfn <- paste(outdir, "Figure1.4_Union_annot.est.png", sep="")
png(figfn, width=520, height=300, res=120)
p
dev.off()




###############################
### summary number of eGene ###
###############################
  
options <- c("torus_dtss", "torus_peak_union")
    
###
for ( option in options){

###    
outdir <- paste("./3_summary.outs/", option, "/summary/", sep="")
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)


condition <- "ALOFT"
fn <- paste("./torus_output/", option,  "/", condition, ".rst", sep="") 
res <- read.table(fn)
    
res <- res%>%dplyr::rename(lfdr=V3)%>%arrange(lfdr)
FDR <- cumsum(res$lfdr)/1:nrow(res)
res$FDR <- FDR
names(res) <- c("order", "gene", "lfdr", "is_egene", "fdr")  
##
opfn <- paste(outdir, condition, ".rst", sep="")
write.table(res, opfn, quote=F, row.names=F, sep="\t")   

##  
}


###
### summarize number of eGenes

summ <- map_dfr(options, function(ii){    
  ###    
  outdir <- paste("./3_summary.outs/", option, "/summary/", sep="")
  fn <- paste(outdir, condition, ".rst", sep="")
  est <- read.table(fn, header=T)
  est2 <- est%>%filter(fdr<0.1)  
  df2 <- data.frame(ntest=length(unique(est$gene)), negene=length(unique(est2$gene)), options=ii)
  df2  
})



