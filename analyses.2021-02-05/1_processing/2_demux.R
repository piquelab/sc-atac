source("../LibraryPackage.R")

### read results of demuxlet
if(FALSE){
basefolder <- "/nfs/rprdata/julong/sc-atac/demux.2021-01-23/demuxOut/"
expNames <- unique(gsub("\\..*", "", list.files(basefolder, "^SCAIP*")))
folders <- paste(basefolder, expNames, ".out.best", sep="")
names(folders) <- expNames

###
demux <- map_dfr(expNames, function(ii){
  fn <- folders[ii]
  dd <- fread(fn)
  dd <- dd%>%mutate(NEW_BARCODE=paste(ii, "_", BARCODE, sep=""), EXP=ii, treats=gsub("SCAIP6[AB]-ATAC-", "", EXP))
  dd
})
### output
opfn <- "./outs/2_demuxlet.rds" ## 115,812 barcodes
write_rds(demux, opfn)
}

### subset SNG cells
if(FALSE){
atac <- read_rds("./outs/1_seurat.merge.rds")
demux <- read_rds("./outs/2_demuxlet.rds")#%>%dplyr::filter(DROPLET.TYPE=="SNG")
meta <- atac@meta.data
meta <- meta%>%left_join(demux, by=c("barcode"="NEW_BARCODE"))
meta0 <- meta%>%dplyr::filter(DROPLET.TYPE=="SNG")
rownames(meta0) <- meta0$barcode
##
atac2 <- subset(atac, cells=meta0$barcode)
atac2@meta.data <- meta0
opfn <- "./outs/2_seurat.merge.SNG.rds"
write_rds(atac2, opfn)
}

#######################
### summary results ###
#######################
### (1)
if(FALSE){
demux <- read_rds("./outs/2_demuxlet.rds") 
demux <- demux%>%mutate(BATCH=gsub("-.*","",EXP))                
demux <- demux%>%dplyr::select("SNG.BEST.GUESS", "EXP", "BATCH", "DROPLET.TYPE")
dd <- demux%>%
   group_by(EXP, DROPLET.TYPE)%>%
   summarise(ncell=n(),.groups="drop")

###         
fig0 <- ggplot(dd)+
   geom_bar(stat="identity", position=position_fill(reverse=F), aes(x=EXP, y=ncell, fill=DROPLET.TYPE))+
   theme_bw()+
   theme(legend.title=element_blank(),
         axis.title=element_blank(), 
         axis.text.x=element_text(angle=90, hjust=1, size=8))   
                
png("./figures/Figure2.1.png", width=600, height=500, res=120)
fig0
dev.off()

### (2)
demux <- read_rds("./outs/2_demuxlet.rds") 
demux <- demux%>%dplyr::filter(DROPLET.TYPE=="SNG")
dd <- demux%>%
   group_by(SNG.BEST.GUESS)%>%
   summarise(ncell=n(),.groups="drop")
fig1 <- ggplot(dd)+
   geom_bar(stat="identity", aes(x=SNG.BEST.GUESS, y=ncell), fill="#1c9099")+
   ggtitle("Number of cells")+
   theme_bw()+
   theme(legend.title=element_blank(),
         axis.title=element_blank(), 
         axis.text.x=element_text(angle=90, hjust=1, size=8),
         plot.title=element_text(hjust=0.5))
png("./figures/Figure2.2_IND.png", width=600, height=500, res=120)
fig1
dev.off()  

### (3)
demux <- read_rds("./outs/2_demuxlet.rds") 
demux <- demux%>%dplyr::filter(DROPLET.TYPE=="SNG")
dd <- demux%>%
   group_by(EXP)%>%
   summarise(ncell=n(),.groups="drop")
fig1 <- ggplot(dd)+
   geom_bar(stat="identity", aes(x=EXP, y=ncell), fill="#1c9099")+
   ggtitle("Number of cells")+
   theme_bw()+
   theme(legend.title=element_blank(),
         axis.title=element_blank(), 
         axis.text.x=element_text(angle=90, hjust=1, size=8),
         plot.title=element_text(hjust=0.5))
png("./figures/Figure2.3_EXP.png", width=500, height=500, res=120)
fig1
dev.off() 

###
demux <- read_rds("./outs/2_demuxlet.rds")
fig1 <- ggplot(demux, aes(x=NUM.SNPS))+
        geom_histogram(fill="grey70", color="grey30", position="identity")+
        xlab("NUM.SNPs")+theme_bw()
 
fig2 <- ggplot(demux, aes(x=NUM.READS))+
        geom_histogram(fill="grey70", color="grey30", position="identity")+
        xlab("NUM.Reads")+theme_bw() 
        
png("./figures/Figure2.4_hist.png", width=600, height=400, res=120)
plot_grid(fig1, fig2,ncol=2)
dev.off()  
}  






