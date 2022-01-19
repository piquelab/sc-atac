###
library(tidyverse)
## library(parallel)
library(data.table)
## library(purrr)
library(GenomicRanges)
library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(Signac)
library(SeuratWrappers)
library(SeuratObject)
## library(cicero, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
## library(monocle3)
###
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(viridis)
theme_set(theme_grey())

outdir <- "2_demux.outs"
if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)

### read results of demuxlet
basefolder <- "/nfs/rprdata/julong/sc-atac/demux.2021-01-23/demuxOut/"
expNames <- unique(gsub("\\..*", "", list.files(basefolder, "^SCAIP*")))
folders <- paste(basefolder, expNames, ".out.best", sep="")
names(folders) <- expNames

###
demux <- map_dfr(expNames, function(ii){
  fn <- folders[ii]
  dd <- fread(fn)
  dd <- dd%>%mutate(BARCODE=gsub("-1","",BARCODE),
      NEW_BARCODE=paste(ii, "_", BARCODE, sep=""),
      EXP=ii, treats=gsub("SCAIP6[AB]-ATAC-", "", EXP))
  dd
})
### output
opfn <- "./2_demux.outs/1_demux.rds" ## 115,812 barcodes
write_rds(demux, opfn)



#########################################
### subset SNG cells of seurat object ###
#########################################
rm(list=ls())
atac <- read_rds("./1_Merge.outs/1_seurat.merge.rds")
atac <- RenameCells(atac, new.names=gsub("-1","",Cells(atac)))
demux <- read_rds("./2_demux.outs/1_demux.rds")%>%
   dplyr::select(NEW_BARCODE, SNG.BEST.GUESS, DROPLET.TYPE)
meta <- atac@meta.data%>%mutate("NEW_BARCODE"=gsub("-1", "", barcode))
meta <- meta%>%left_join(demux, by="NEW_BARCODE")
rownames(meta) <- meta$NEW_BARCODE
atac <- AddMetaData(atac, meta)

##
atac2 <- subset(atac, subset=DROPLET.TYPE=="SNG")
opfn <- "./2_demux.outs/2_seurat.merge.SNG.rds"
write_rds(atac2, opfn)



#######################
### summary results ###
#######################

demux <- read_rds("./2_demux.outs/1_demux.rds") 
demux <- demux%>%mutate(BATCH=gsub("-.*","",EXP))                
## demux <- demux%>%dplyr::select("SNG.BEST.GUESS", "EXP", "BATCH", "DROPLET.TYPE")

##(1)
dd <- demux%>%
   group_by(EXP, DROPLET.TYPE)%>%
   summarise(ncell=n(),.groups="drop")        
fig1 <- ggplot(dd)+
   geom_bar(stat="identity",
            position=position_fill(reverse=F),
            aes(x=EXP, y=ncell, fill=DROPLET.TYPE))+
   theme_bw()+
   theme(legend.title=element_blank(),
         axis.title=element_blank(), 
         axis.text.x=element_text(angle=90, hjust=1, size=8))   
                
png("./2_demux.outs/Figure1.1_droplet.png", width=600, height=500, res=120)
fig1
dev.off()


### (2)
demux2 <- demux%>%dplyr::filter(DROPLET.TYPE=="SNG")
dd <- demux2%>%
   group_by(SNG.BEST.GUESS)%>%
   summarise(ncell=n(),.groups="drop")
fig2 <- ggplot(dd)+
   geom_bar(stat="identity", aes(x=SNG.BEST.GUESS, y=ncell), fill="#1c9099")+
   ggtitle("Number of cells")+
   theme_bw()+
   theme(legend.title=element_blank(),
         axis.title=element_blank(), 
         axis.text.x=element_text(angle=90, hjust=1, size=8),
         plot.title=element_text(hjust=0.5))
png("./2_demux.outs/Figure1.2_IND.png", width=600, height=500, res=120)
fig2
dev.off()


### (3)
dd <- demux2%>%
   group_by(EXP)%>%
   summarise(ncell=n(),.groups="drop")
fig3 <- ggplot(dd)+
   geom_bar(stat="identity", aes(x=EXP, y=ncell), fill="#1c9099")+
   ggtitle("Number of cells")+
   theme_bw()+
   theme(legend.title=element_blank(),
         axis.title=element_blank(), 
         axis.text.x=element_text(angle=90, hjust=1, size=8),
         plot.title=element_text(hjust=0.5))
png("./2_demux.outs/Figure1.3_EXP.png", width=500, height=500, res=120)
fig3
dev.off() 


### (4)
fig4.1 <- ggplot(demux, aes(x=NUM.SNPS))+
        geom_histogram(fill="grey70", color="grey30", position="identity")+
        xlab("NUM.SNPs")+theme_bw()
 
fig4.2 <- ggplot(demux, aes(x=NUM.READS))+
        geom_histogram(fill="grey70", color="grey30", position="identity")+
        xlab("NUM.Reads")+theme_bw() 
        
png("./2_demux.outs/Figure1.4_hist.png", width=600, height=400, res=120)
plot_grid(fig4.1, fig4.2, ncol=2)
dev.off()  



######################
### summary reads  ###
#####################

fn <- "./2_demux.outs/2_seurat.merge.SNG.rds"
atac<- read_rds(fn)

## counts <- GetAssayData(sc, slot="counts")
## nCount_ATAC <- colSums(counts)
## nFeature_ATAC <- colSums(counts>0)
###
meta <- atac@meta.data

meta2 <- meta%>%mutate(treats=gsub(".*-ATAC-|_.*", "", NEW_BARCODE),
                       treat2=gsub("-EtOH", "", treats),
                       EXP=gsub("_.*","",NEW_BARCODE))


dd2 <- meta2%>%group_by(SNG.BEST.GUESS, treat2)%>%
   summarise(ncell=n(), reads=mean(nCount_ATAC), ngene=mean(nFeature_ATAC),.groups="drop")



xx2 <- meta2%>%group_by(EXP)%>%
   summarise(ncell=n(), reads=mean(nCount_ATAC), ngene=mean(nFeature_ATAC),.groups="drop")

tmp <- dd2%>%
           group_by(treat2)%>%
           summarise(nind=n(), ncell=median(ncell), reads=median(reads),ngene=median(ngene),.groups="drop")


col1 <- c("CTRL"="#828282",
   "LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
   "PHA"="#a6cee3", "PHA-DEX"="#1f78b4")
lab1 <- c("CTRL"="CTRL",
   "LPS"="LPS", "LPS-DEX"="LPS+DEX",
   "PHA"="PHA", "PHA-DEX"="PHA+DEX")


fig1 <- ggplot(dd2, aes(x=treat2, y=ncell, fill=treat2))+
   geom_violin()+xlab("")+ylab("")+
   ggtitle("#Cells per individual")+
   scale_fill_manual(values=col1)+
   scale_x_discrete(labels=lab1)+
   theme_bw()+
   theme(legend.position="none",
   plot.title=element_text(hjust=0.5),
   axis.text.x=element_text(angle=-90, hjust=0, vjust=0.5))

###
fig2 <- ggplot(dd2,aes(x=treat2, y=reads, fill=treat2))+
   geom_violin()+xlab("")+ylab("")+
   ggtitle("#UMIs per cell")+
   scale_fill_manual(values=col1)+
   scale_x_discrete(labels=lab1)+
   theme_bw()+
   theme(legend.position="none",
         plot.title=element_text(hjust=0.5),
         axis.text.x=element_text(angle=-90, hjust=0, vjust=0.5))

###
fig3 <- ggplot(dd2,aes(x=treat2, y=ngene, fill=treat2))+
    geom_violin()+xlab("")+ylab("")+
    ggtitle("#Genes per cell")+
    scale_fill_manual(values=col1)+
    scale_x_discrete(labels=lab1)+
    theme_bw()+
    theme(legend.position="none",
          plot.title=element_text(hjust=0.5),
          axis.text.x=element_text(angle=-90, hjust=0, vjust=0.5))

## png("./2_kb2_output/Figure2.5_violin.png", width=800, height=500, res=120)
pdf("./2_demux.outs/Figure2.1_violin.pdf", width=8, height=5)
print(plot_grid(fig1, fig2, fig3, labels="AUTO", label_fontface="plain", label_x=0.1,  ncol=3))
dev.off()

