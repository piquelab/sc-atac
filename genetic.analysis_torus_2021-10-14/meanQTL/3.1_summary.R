##
library(tidyverse)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(viridis)


#################################################################
### forest plot for cell-type and response motids annotation  ###
#################################################################

#### final version

rm(list=ls())

pct0 <- 0.02

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
      scale_x_continuous("log odds ratio", breaks=seq(-15, 8, 3), limits=c(-15,8)) +
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
