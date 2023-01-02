##
library(tidyverse)
library(data.table)
library(qvalue)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(RColorBrewer)
library(viridis)

rm(list=ls())
 
outdir <- "./3_summary.outs/"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)


################################################################################################
### 1. combine annotation, peaks, cell-type motifs, treatment motifs in different cell-types ###
################################################################################################

fn <- "./torus_output/Combine/Asthma_combine2_all.est"
est <- read.table(fn)
est2 <- est[2:22,]
est2 <- est2%>%mutate(order=as.numeric(1:nrow(est2)),
   category=as.character(gsub(".1", "", V1)),
   gr=ifelse(grepl("_", category), gsub(".*_", "", category), category))%>%
   dplyr::rename("odds"="V2", "CI_lower"="V3", "CI_upper"="V4")
est2 <- est2%>%mutate(category2=fct_reorder(category, order))


col2 <- c("Bcell"="#4daf4a", "Monocyte"="#984ea3","NKcell"="#aa4b56", "Tcell"="#ffaa00",
  "LPS"="#fb9a99", "LPS-DEX"="#e31a1c", "PHA"="#a6cee3", "PHA-DEX"="#1f78b4", "peaking"="#e7298a")

p1 <- ggplot(est2, aes(x=odds, y=category2, color=factor(gr)))+
    geom_errorbarh(aes(xmax=CI_upper, xmin=CI_lower), size=0.5, height=0.2)+
    geom_point(shape=19, size=0.5)+
    geom_vline(aes(xintercept=0), size=0.25, linetype="dashed")+
    ##scale_y_discrete(labels=ylab)+
    scale_x_continuous("log odds ratio", breaks=seq(-4, 6, by=2), limits=c(-5, 7))+
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
 
figfn <- paste(outdir, "Figure1.1_combine2_annot.est.png", sep="")
png(figfn, width=520, height=500, res=120)
p1
dev.off()


######################################################
### 2. peaks, cell-type motifs and response motifs ###
######################################################

fn <- "./torus_output/Combine/Asthma_Union_all.est"
est <- read.table(fn)
est2 <- est[2:4,]
est2 <- est2%>%mutate(order=as.numeric(1:nrow(est2)))%>%
   dplyr::rename("peaking_d"="V1", "odds"="V2", "CI_lower"="V3", "CI_upper"="V4")

p2 <- ggplot(est2, aes(x=odds, y=as.character(order), color=factor(order)))+
    geom_errorbarh(aes(xmax=CI_upper, xmin=CI_lower), size=0.5, height=0.2)+
    geom_point(shape=19, size=0.5)+
    geom_vline(aes(xintercept=0), size=0.25, linetype="dashed")+
    scale_y_discrete(labels=c("1"="Peaks", "2"="Cell-type motifs", "3"="Treatment motifs"))+
    scale_x_continuous("log odds ratio", breaks=seq(-2,6,by=2), limits=c(-2, 7))+
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
 
figfn <- paste(outdir, "Figure1.2_union_annot.est.png", sep="")
png(figfn, width=480, height=300, res=120)
p2
dev.off()



## ##########################
## ### meta T-cells scaip ###
## ##########################
## fn <- "/nfs/rprdata/julong/sc-atac/genetic.analysis_torus_2021-10-14/eQTL/datasets.txt" 
## datasets <- read.table(fn)$V1
## est3 <- map_dfr(16:20, function(i){
##     ##
##     condition <- datasets[i]
##     fn <- paste("/nfs/rprdata/julong/sc-atac/genetic.analysis_torus_2021-10-14/eQTL/torus_output/torus_peak_union/", condition, ".est", sep="")
##     est <- read.table(fn)
##     est <- est[2:4,]
##     est <- est%>%mutate(conditions=condition,
##        MCls=gsub("_.*", "", conditions),
##        order=gsub("peaking.", "", V1))
## ###
##     est
## })
## ###
## est3 <- est3%>%dplyr::rename("odds"="V2", "CI_lower"="V3", "CI_upper"="V4")%>%
##     mutate(se=abs(odds-CI_lower)/1.96)

## DF <- NULL
## for (i in 1:3){
##    ##
##    est0 <- est3%>%dplyr::filter(order==i)    
##    b <- est0$odds 
##    v <- (est0$se)^2
##    ##
##    bmeta <- sum((1/v)*(1/sum(1/v))*b)
##    vmeta <- 1/sum(1/v)
##    se <- sqrt(vmeta) 
##    df2 <- data.frame(order=i, b=bmeta, se=sqrt(vmeta), CI_lower=bmeta-1.96*se, CI_upper=bmeta+1.96*se)
##    DF <- rbind(DF, df2)
## }

## DF$MCls <- "Tcell"
## DF <- DF%>%dplyr::select(odds=b, CI_lower, CI_upper, order, MCls)

## ####
## #### ALOFT
## fn <- "./torus_output/torus_peak_combineNew/ALOFT_Union.est"
## est <- read.table(fn)
## est2 <- est[2:4,]
## est2 <- est2%>%mutate(order=as.numeric(1:nrow(est2)), MCls="Bulk")%>%
##    dplyr::rename("peaking_d"="V1", "odds"="V2", "CI_lower"="V3", "CI_upper"="V4")%>%
##    dplyr::select(odds, CI_lower, CI_upper, order, MCls)

## plotDF <- rbind(DF, est2)%>%
##     mutate(comb=paste(order, MCls, sep="_"))

## ###
## ###
## p <- ggplot(plotDF,
##     aes(x=odds, y=as.character(order), color=factor(order), group=factor(MCls), alpha=factor(MCls)))+
##     geom_errorbarh(aes(xmax=CI_upper, xmin=CI_lower), size=0.5, height=0.2, position=position_dodge(0.5))+
##     geom_point(shape=19, size=0.5, position=position_dodge(0.5))+
##     geom_vline(aes(xintercept=0), size=0.25, linetype="dashed")+
##     scale_y_discrete(labels=c("1"="Peaks", "2"="Cell-type motifs", "3"="Treatment motifs"))+
##     scale_x_continuous("log odds ratio", limits=c(-1, 4.5))+
##     scale_colour_manual("", values=c("1"="#7570b3", "2"="#1b9e77", "3"="#e7298a"),
##        labels=c("1"="Peaks", "2"="Cell-type motifs", "3"="Treatment motifs"), guide="none")+
##    scale_alpha_manual("", values=c("Tcell"=1, "Bulk"=0.5),
##        labels=c("Tcell"="T-cells (SCAIP)", "Bulk"="Whole blood (ALOFT)"),               
##        guide=guide_legend(override.aes=list(color="#e7298a")))+ 
##     theme_bw()+
##     theme(legend.title=element_blank(),
##           legend.key.size=grid::unit(1, "lines"),
##           axis.title.y=element_blank())
##           ##axis.title.x=element_text(size=10),
##           ###strip.text=element_text(size=12),
##           ##plot.title=element_text(hjust=0.5, size=12),
##           ## axis.text.x=element_text(size=9, angle=-90, hjust=0, vjust=0.5),
##           ###axis.text=element_text(size=9))
 
## figfn <- paste(outdir, "Figure1.4_Union_annot.est.png", sep="")
## png(figfn, width=520, height=300, res=120)
## p
## dev.off()






