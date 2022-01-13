##
library(tidyverse)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(viridis)


###############
### summary ###
###############

###
outdir <- "./torus_output/torus_test/"

### table of egenes
conditions <- read.table("datasets.txt", header=F)$V1
summ <- map_dfr(conditions, function(ii){
   ###
   fn <- paste(outdir, ii, ".egene.rst", sep="")
   res <- read.table(fn, header=F)
   ###
   df <- data.frame(condition=ii, ngene=sum(res$V4))
   df                 
})

opfn <- paste(outdir, "zzz_summary.csv", sep="")
write.csv(summ, opfn, row.names=F)

### eGene list
egene <- map_dfr(conditions, function(ii){
   ###
   fn <- paste(outdir, ii, ".egene.rst", sep="")
   res <- read.table(fn, header=F)
   ###
   res2 <- res%>%dplyr::filter(V4==1)
})
###
opfn <- paste(outdir, "zzz_egene.csv", sep="")
write.csv(egene, opfn, row.names=F)

length(unique(egene$V2))



##################
### forestplot ###
##################

outdir <- "./3_summary.outs/"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)


datasets <- read.table("datasets.txt")

est2 <- map_dfr(1:20, function(i){
###
  condition <- datasets$V1[i]
  fn <- paste("./torus_output/torus_tss/", condition, ".est", sep="") 
  est <- read.table(fn)
  est <- est%>%mutate(order=1:nrow(est), conditions=condition,
      MCls=gsub("_.*", "", conditions),
      treats=gsub("-EtOH", "", gsub(".*_", "", conditions)),
      label.name=gsub("]", "", gsub(".*\\[", "", V1)),
      dtss=as.numeric(gsub("kb","",label.name)))
})


###
tmp <- est2%>%filter(conditions=="Bcell_CTRL")
ylabel.names <- as.character(tmp$label.name)
names(ylabel.names) <- tmp$order

p <- ggplot(est2, aes(x=V2, y=as.factor(order),color=factor(treats)))+
      geom_errorbarh(aes(xmax=V4, xmin=V3), size=0.5, height=0.2)+
      geom_point(shape=19, size=0.5)+
      geom_vline(aes(xintercept=0), size=0.25, linetype="dashed")+
      xlab("log odds ratio")+
      scale_y_discrete(labels=ylabel.names)+
      scale_colour_manual("Treatment",
         values=c("CTRL"="#828282", "LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
                  "PHA"="#a6cee3", "PHA-DEX"="#1f78b4"),
         guide=guide_legend(override.aes=list(size=1)) )+
     
      facet_wrap(~MCls, ncol=4)+
      theme_bw()+
      theme(## plot.title=element_text(hjust=0.5, size=10),
            axis.title.x=element_text(size=10),
            axis.title.y=element_blank(),
            axis.text=element_text(size=8))
###
figfn <- paste(outdir, "torus_tss/Figure0_tss", ".est.png", sep="")
png(figfn, width=800, height=500, res=120)
print(p)
dev.off()



###
MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
for (ii in MCls){
cat(ii,"\n")    
###    
p <- ggplot(est2%>%filter(MCls==ii), aes(x=V2, y=as.factor(order),color=factor(treats)))+
      geom_errorbarh(aes(xmax=V4, xmin=V3), size=0.5, height=0.2)+
      geom_point(shape=19, size=0.5)+
      geom_vline(aes(xintercept=0), size=0.25, linetype="dashed")+
      xlab("log odds ratio")+ylab("Distance to TSS")+
      scale_y_discrete(labels=ylabel.names)+
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
figfn <- paste(outdir, "torus_tss/Figure0.1_", ii, ".tss.est.png", sep="")
png(figfn, width=800, height=500, res=120)
print(p)
dev.off()
}###




########################################
#### forest plot for response motifs ###
########################################

outdir <- "./3_summary.outs/torus_motif2/"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)


datasets <- read.table("datasets.txt")
response <- read.table("datasets_motif.txt")$V1
###
est2 <- map_dfr(1:nrow(datasets), function(i){
###
  condition <- datasets$V1[i]
  est <- map_dfr(1:length(response), function(j){
      contrast <- response[j]
      fn <- paste("./torus_output/torus_motif2/", contrast, "/",  condition, ".est", sep="") 
      est0 <- read.table(fn)
      est0[2,]
  })
  est <- est%>%mutate(conditions=condition,
      MCls=gsub("_.*", "", conditions), treats=gsub(".*_|-EtOH", "", condition),
      motif_conditions=gsub("\\..*", "", V1),
      motif_MCls=gsub("_.*", "", motif_conditions),
      motif_contrasts=gsub(".*_", "", motif_conditions))
    
###
})

####
MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
for (ii in MCls){
cat(ii,"\n")    
###    
p <- ggplot(est2%>%filter(MCls==ii), aes(x=V2, y=motif_conditions, color=factor(motif_MCls)))+
      geom_errorbarh(aes(xmax=V4, xmin=V3), size=0.5, height=0.2)+
      geom_point(shape=19, size=0.5)+
      geom_vline(aes(xintercept=0), size=0.25, linetype="dashed")+
      xlab("log odds ratio")+
    ## scale_y_discrete(labels=ylabel.names)+
      scale_colour_manual(
          values=c("Bcell"="#4daf4a", "Monocyte"="#984ea3",
                   "NKcell"="#aa4b56", "Tcell"="#ffaa00"))+
      ggtitle(ii)+
      facet_wrap(~treats, ncol=5)+
      theme_bw()+
      theme(plot.title=element_text(hjust=0.5, size=10),
            legend.position="none",
            axis.title.x=element_text(size=10),
            axis.title.y=element_blank(),
            axis.text=element_text(size=8))
###
figfn <- paste(outdir, "Figure0.1_", ii, ".tss.est.png", sep="")
png(figfn, width=800, height=500, res=120)
print(p)
dev.off()
}###



######################
### positive motif ###
######################

outdir <- "./3_summary.outs/torus_motif3/"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)


datasets <- read.table("datasets.txt")
response <- read.table("datasets_motif.txt")$V1
###
est2 <- map_dfr(1:nrow(datasets), function(i){
###
  condition <- datasets$V1[i]
  est <- map_dfr(1:length(response), function(j){
      contrast <- response[j]
      fn <- paste("./torus_output/torus_motif3/", contrast, "/",  condition, ".est", sep="") 
      est0 <- read.table(fn)
      est0[2,]
  })
  est <- est%>%mutate(conditions=condition,
      MCls=gsub("_.*", "", conditions), treats=gsub(".*_|-EtOH", "", condition),
      motif_conditions=gsub("\\..*", "", V1),
      motif_MCls=gsub("_.*", "", motif_conditions),
      motif_contrasts=gsub(".*_", "", motif_conditions))
    
###
})

####
MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
for (ii in MCls){
cat(ii,"\n")    
###    
p <- ggplot(est2%>%filter(MCls==ii), aes(x=V2, y=motif_conditions, color=factor(motif_MCls)))+
      geom_errorbarh(aes(xmax=V4, xmin=V3), size=0.5, height=0.2)+
      geom_point(shape=19, size=0.5)+
      geom_vline(aes(xintercept=0), size=0.25, linetype="dashed")+
      xlab("log odds ratio")+
    ## scale_y_discrete(labels=ylabel.names)+
      scale_colour_manual(
          values=c("Bcell"="#4daf4a", "Monocyte"="#984ea3",
                   "NKcell"="#aa4b56", "Tcell"="#ffaa00"))+
      ggtitle(ii)+
      facet_wrap(~treats, ncol=5)+
      theme_bw()+
      theme(plot.title=element_text(hjust=0.5, size=10),
            legend.position="none",
            axis.title.x=element_text(size=10),
            axis.title.y=element_blank(),
            axis.text=element_text(size=8))
###
figfn <- paste(outdir, "Figure0.1_", ii, ".tss.est.png", sep="")
png(figfn, width=800, height=500, res=120)
print(p)
dev.off()
}###







## for (i in 1:20){
## ###
##    condition <- datasets$V1[i]
##    fn <- paste("./torus_output/torus_tss/", condition, ".est", sep="") 
##    est <- read.table(fn)## %>%mutate(V2=exp(V2),V3=exp(V3),V4=exp(V4))
##    est$order <- 1:nrow(est)

##    est2 <- est

##    ### 
##    ylab.name <- gsub("]", "", gsub(".*\\[", "", est2$V1))   
##    names(ylab.name) <- as.character(est2$order)
##    fig <- ggplot(est2, aes(x=V2, y=as.factor(order)))+
##       geom_errorbarh(aes(xmax=V4, xmin=V3), size=0.5, height=0.2)+
##       geom_point(shape=19, size=2.5)+
##       geom_vline(aes(xintercept=0), size=0.25, linetype="dashed")+
##       xlab("log odds ratio")+
##       scale_y_discrete(labels=ylab.name)+ 
##       ggtitle(condition)+
##       theme_bw()+
##       theme(plot.title=element_text(hjust=0.5, size=10),
##             axis.title.x=element_text(size=10),
##             axis.title.y=element_blank(),
##             axis.text.x=element_text(size=8))
   
##    figfn <- paste(outdir, "Figure", i, "_", condition, ".est.png", sep="")
##    png(figfn, width=400, height=600, res=120)
##    print(fig)
##    dev.off()
## }
