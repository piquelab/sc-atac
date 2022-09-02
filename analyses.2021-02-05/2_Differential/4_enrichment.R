
##
library(tidyverse)
library(annotables)
library(clusterProfiler)
library(org.Hs.eg.db)
library(EnsDb.Hsapiens.v75)
library(ChIPseeker,lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")

##
library(ggplot2)
library(cowplot)
library(grid)
library(gridExtra)
library(ggExtra)
library(RColorBrewer)
library(gtable)
library(ggsignif)
library(pheatmap)
library(corrplot)
library(viridis)

## theme_set(theme_grey())

rm(list=ls())

###
outdir <- "./4_enrichment.outs/"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)


###
### select peaks at at least 2% cells for each cell-type
atac <- read_rds("../1_processing/5.1_reCallPeak.outs/3_scATAC.annot.rds")
MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell", "DC")
peakSel <- lapply(MCls, function(ii){
    ##
    atac2 <- subset(atac, subset=MCls==ii)
    count <- atac2@assays$ATAC@counts
    rpz <- rowMeans(count>0)
    peakSel2 <- rownames(count)[rpz>0.02]
    peakSel2
})
peakSel <- do.call(c, peakSel)
peakSel <- unique(peakSel)


###
###

###
res <- read_rds("./1.3_DiffPeak.outs/3.0_DESeq_indi.results.rds")
## res <- res%>%filter(gene%in%peakSel)

peakAnno <- read_rds("2.2_compareRNAandATAC.outs/2_annot.ChIPseeker.rds")%>%
   as.data.frame()%>%
   mutate(chr=gsub("chr", "", seqnames),
          peak_region=paste(chr,start,end,sep="-"),
          dis=abs(distanceToTSS))%>%
   dplyr::select(peak_region, geneId, distanceToTSS, dis)# flank_geneIds) 
    

###
### background gene
#gene <- unique(unlist(str_split(peakAnno$flank_geneIds, ";")))
gene <- unique(peakAnno$geneId)
x <- grch37%>%dplyr::filter(ensgene%in%gene, grepl("protein_coding", biotype))
geneBG <- bitr(x$ensgene, fromType="ENSEMBL", toType=c("ENTREZID","SYMBOL"),
               OrgDb=org.Hs.eg.db)

###
### gene cluster for enrichment data analysis
res2 <- res%>%as.data.frame()%>%dplyr::rename("peak_region"="gene")%>%
    mutate(comb=paste(MCls, contrast, sep="_"), direction=ifelse(estimate>0, "Up", "Down"))%>%
    dplyr::filter(abs(estimate)>0.5, p.adjusted<0.1, MCls!="DC")
res2 <- res2%>%left_join(peakAnno, by="peak_region")

##
## comb <- sort(unique(res2$comb))
## df0 <- map_dfr(comb, function(ii){
##    tmp <- res2%>%dplyr::filter(comb==ii)%>%mutate(dis=abs(distanceToTSS))%>%
##       group_by(geneId)%>%top_n(-1, dis)
##       dplyr::filter(dis==min(dis,na.rm=T))%>%ungroup()
##    tmp
## })

df0 <- res2%>%group_by(contrast, MCls, geneId)%>%top_n(-1, dis)%>%ungroup()%>%as.data.frame()

df2 <- df0%>%dplyr::filter(dis<100000)

x <- bitr(df2$geneId, fromType="ENSEMBL",
          toType=c("ENTREZID","SYMBOL"), OrgDb=org.Hs.eg.db)

geneCluster <- df2%>%inner_join(x, by=c("geneId"="ENSEMBL"))


###
### GO enrichment
# "GO analysis for direction (up and down) directly", "\n")
cg <- compareCluster(ENTREZID~contrast+MCls+direction,
    data=geneCluster,
    universe=geneBG$ENTREZID,
    fun="enrichGO",
    OrgDb="org.Hs.eg.db",
    pvalueCutoff=1,
    qvalueCutoff=1,
    ont="ALL",
    minGSSize=0,
    maxGSSize=1000)

#cg0 <- as.data.frame(cg)
cg <- cg%>%mutate(Cluster1=gsub("\\.((Down)|(Up))", "", Cluster))
write_rds(cg,"./4_enrichment.outs/1_enrichGO.rds")


####################
### show figures ###
####################

library(ggtext, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
library(glue)

cg <- read_rds("./4_enrichment.outs/1_enrichGO.rds")

cl <- paste( rep(c("LPS", "PHA", "LPS+DEX", "PHA+DEX"),each=4),
   rep(c("Bcell", "Monocyte", "NKcell", "Tcell"), times=4), sep=".")
cl <- paste(rep(cl,times=2), rep(c("Up","Down"), each=16), sep=".")

cluster2 <- 1:32
names(cluster2) <- cl


## ii <- paste(rep(c("A","B","C","D"),each=4),rep(1:4,times=4), sep="")
## cl2 <- paste(rep(c("X","Y"),each=16), rep(ii,times=2), sep=".")
## cluster2 <- setNames(cl2, cl)
## lab2 <- setNames(gsub("-","+",cl),cl2)
col1 <- c("LPS"="#fb9a99", "LPS+DEX"="#e31a1c", "PHA"="#a6cee3", "PHA+DEX"="#1f78b4")
col2 <- c("Bcell"="#4daf4a", "Monocyte"="#984ea3", "NKcell"="#aa4b56", "Tcell"="#ffaa00")


cg <- cg%>%
  mutate(Cluster=gsub("-", "+", Cluster),
         contrast=gsub("-", "+", contrast),
         col.contrast=col1[contrast],
         col.MCls=col2[MCls],
         ClusterValue=as.numeric(cluster2[Cluster]),
         ClusterNew=glue("<i style='color:{col.contrast}'>{contrast}.<i style='color:{col.MCls}'>{MCls}.<i style='color:black'>{direction}"),
         ClusterNew=fct_reorder(ClusterNew, ClusterValue),
         maxGSSize=as.numeric(gsub("/.*", "", BgRatio)),
         ngene=as.numeric(gsub(".*/", "", GeneRatio)) )

cg2 <- cg%>%dplyr::filter(maxGSSize>10, maxGSSize<500, ngene>10, p.adjust<0.1)

fig1 <- enrichplot::dotplot(cg2, x="ClusterNew", showCategory=5)+
   theme(axis.title=element_blank(),
         axis.text.x=element_markdown(angle=45, hjust=1,size=16),
         axis.text.y=element_text(size=10))

### 
figfn <- "./4_enrichment.outs/Figure1.1_GO.png"
png(figfn, width=3500, height=2500, res=160)
print(fig1)
dev.off()


###
figfn <- "./4_enrichment.outs/Figure1.1_GO.pdf"
pdf(figfn, width=18, height=15)
print(fig1)
dev.off()


#####################
### KEGG analysis ###
#####################

ck <- compareCluster(ENTREZID~contrast+MCls+direction,
    data=geneCluster,
    universe=geneBG$ENTREZID,
    fun="enrichKEGG",
    ## OrgDb="org.Hs.eg.db",
    pvalueCutoff=1,
    qvalueCutoff=1,
    ## ont="ALL",
    minGSSize=0,
    maxGSSize=1000)

ck <- ck%>%mutate(Cluster1=gsub("\\.((Down)|(Up))", "", Cluster))
write_rds(ck,"./4_enrichment.outs/2_enrichKEGG.rds")


###
### figures

ck <- read_rds("./4_enrichment.outs/2_enrichKEGG.rds")


cl <- paste( rep(c("LPS", "PHA", "LPS+DEX", "PHA+DEX"),each=4),
   rep(c("Bcell", "Monocyte", "NKcell", "Tcell"), times=4), sep=".")
cl <- paste(rep(cl,times=2), rep(c("Up","Down"), each=16), sep=".")

cluster2 <- 1:32
names(cluster2) <- cl


## ii <- paste(rep(c("A","B","C","D"),each=4),rep(1:4,times=4), sep="")
## cl2 <- paste(rep(c("X","Y"),each=16), rep(ii,times=2), sep=".")
## cluster2 <- setNames(cl2, cl)
## lab2 <- setNames(gsub("-","+",cl),cl2)
col1 <- c("LPS"="#fb9a99", "LPS+DEX"="#e31a1c", "PHA"="#a6cee3", "PHA+DEX"="#1f78b4")
col2 <- c("Bcell"="#4daf4a", "Monocyte"="#984ea3", "NKcell"="#aa4b56", "Tcell"="#ffaa00")


ck <- ck%>%
  mutate(Cluster=gsub("-", "+", Cluster),
         contrast=gsub("-", "+", contrast),
         col.contrast=col1[contrast],
         col.MCls=col2[MCls],
         ClusterValue=as.numeric(cluster2[Cluster]),
         ClusterNew=glue("<i style='color:{col.contrast}'>{contrast}.<i style='color:{col.MCls}'>{MCls}.<i style='color:black'>{direction}"),
         ClusterNew=fct_reorder(ClusterNew, ClusterValue),
         maxGSSize=as.numeric(gsub("/.*", "", BgRatio)),
         ngene=as.numeric(gsub(".*/", "", GeneRatio)) )

ck2 <- ck%>%dplyr::filter(maxGSSize>10, maxGSSize<500, ngene>10, p.adjust<0.1)

fig2 <- enrichplot::dotplot(ck2, x="ClusterNew", showCategory=5)+
   theme(axis.title=element_blank(),
         axis.text.x=element_markdown(angle=45, hjust=1,size=16),
         axis.text.y=element_text(size=13))

###
figfn <- "./4_enrichment.outs/Figure2.1_KEGG.png"
png(figfn, width=2000, height=1500, res=120)
print(fig2)
dev.off()

##
figfn <- "./4_enrichment.outs/Figure2.1_KEGG.pdf"
pdf(figfn, width=18, height=15)
print(fig2)
dev.off()



#################
### examples  ###
#################

###
odds.fun <-  function(df){
###    
   res <- map_dfr(1:nrow(df), function(i){
      Diff <- as.numeric(df[i, c("Diff.in", "Diff.not")])
      Bg <- as.numeric(df[i, c("Bg.in", "Bg.not")])
      dat <- data.frame(Diff=Diff, Bg=Bg)
      rownames(dat) <- c("in.category", "not.category")
      fish <- fisher.test(dat)
      res0 <- data.frame(odds=as.numeric(fish$estimate),
                         down=fish$conf.int[1],
                         up=fish$conf.int[2])
      res0
   })
###
  df$odds <- res$odds
  df$down <- res$down
  df$up <- res$up  
  df  
}


ExampleGOplot <- function(cg, nbreak){ 

### prepare data    
   ## x <- str_split(cg$GeneRatio, "/", simplify=T)
   ## GeneRatio <- as.numeric(x[,1])/as.numeric(x[,2])
   Drt2 <- c("Up"=1, "Down"=2) 
   cg <- cg%>%mutate(Direction2=Drt2[direction],
      contrast2=paste(Direction2, contrast, sep="."))%>%
      mutate(contrast2=gsub("-", "+", contrast2)) 
   ## cg$size <- rep(1,nrow(cg))
   ## cg$size[GeneRatio>=0.05&GeneRatio<0.15] <- 2
   ## cg$size[GeneRatio>=0.15] <- 3 
   #
   ## cg <- cg%>%drop_na(odds)
   fig0 <- ggplot(cg, aes(x=contrast2, y=MCls))+
      geom_point(aes(size=odds, colour=factor(p2)))+
      scale_x_discrete(labels=c("1.LPS"="LPS.Up", "2.LPS"="LPS.Down",
         "1.LPS+DEX"="LPS+DEX.Up", "2.LPS+DEX"="LPS+DEX.Down",
         "1.PHA"="PHA.Up", "2.PHA"="PHA.Down",
         "1.PHA+DEX"="PHA+DEX.Up", "2.PHA+DEX"="PHA+DEX.Down"))+
      scale_colour_manual(name="p.adjust",
          values=c("1"="grey80", "2"="red"),
          labels=c("1"=">0.1", "2"="<0.1"),
          guide=guide_legend(override.aes=list(size=2), order=1))+ 
      ## scale_colour_gradient(name="p.adjust",                           
      ##    low="blue", high="red", na.value=NA, trans="reverse", n.breaks=5,
      ##    guide=guide_colourbar(order=1))+    #"#ffa500
      scale_size_binned("odds ratio",
         breaks=waiver(), n.breaks=nbreak,                
         guide=guide_bins(show.limits=TRUE, axis=TRUE,
           axis.show=arrow(length=unit(1.5,"mm"), ends="both"),
           keywidth=grid::unit(0.4, "lines"),
           keyheight=grid::unit(0.4, "lines"), order=2))+
      theme_bw()+
      theme(axis.title=element_blank(),
         axis.text.y=element_text(size=9), 
         axis.text.x=element_text(angle=-90, size=9, hjust=0, vjust=0.5), 
         legend.background=element_blank(),
         legend.title=element_text(size=8),
         legend.text=element_text(size=6),
         legend.key.size=grid::unit(0.4, "lines"))
   fig0
}



###

contrast <- c("LPS", "LPS-DEX", "PHA", "PHA-DEX")
MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
rn <- paste(rep(contrast, each=8), rep(rep(MCls, each=2), times=4),
            rep(rep(c("Down", "Up"),times=4), times=4), sep=".")
tmp <- data.frame(contrast=rep(contrast, each=8),
   MCls=rep(rep(MCls, each=2), times=4),
   direction=rep(rep(c("Down", "Up"),times=4), times=4))%>%
   mutate(rn=paste(contrast, MCls, direction, sep="."))

###
###
cg <- read_rds("./4_enrichment.outs/1_enrichGO.rds")




###
### response to LPS

## mycolor <- colorRampPalette(c("blue", "red"))(3)

cg2 <- cg%>%dplyr::filter(Description=="response to lipopolysaccharide")

cg2 <- cg2%>%as.data.frame()%>%
    mutate(Diff.in=as.numeric(gsub("/.*","",GeneRatio)),
           Diff.total=as.numeric(gsub(".*/","",GeneRatio)),
           Diff.not=Diff.total-Diff.in)
cg2 <- cg2%>%
    mutate(Bg.in=as.numeric(gsub("/.*","", BgRatio)),
           Bg.total=as.numeric(gsub(".*/","", BgRatio)),
           Bg.not=Bg.total-Bg.in)

cg2 <- odds.fun(cg2)
cg2 <- cg2%>%dplyr::select(Cluster, p.adjust, odds)%>%full_join(tmp, by=c("Cluster"="rn"))
##
##
## pp <- cg2$p.adjust
## p2 <- rep(NA, length(pp))
## p2[pp<=0.3] <- 1
## p2[pp<=0.1] <- 2
## p2[pp<0.01] <- 3
## cg2$p2 <- as.character(p2)
cg2 <- cg2%>%mutate(p2=ifelse(p.adjust<0.1, 2, 1))


fig <- ExampleGOplot(cg2, nbreak=3)+
    ## scale_colour_manual(name="p.adjust",
    ##    values=c("1"=mycolor[1], "2"=mycolor[2], "3"=mycolor[3]),
    ##    labels=c("1"="0.1~0.2", "2"="0.01~0.1", "3"="~0.01"), na.value=NA,
    ##    guide=guide_legend(override.aes=list(size=2), order=1))+
   ggtitle("response to LPS")+
   theme(plot.title=element_text(hjust=0.5, size=14))

figfn <-"./4_enrichment.outs/Figure3.1_response_to_LPS.png"
png(figfn, width=550, height=400, res=120)
print(fig)
dev.off()


###
### type I interferon
cg2 <- cg%>%dplyr::filter(Description=="type I interferon signaling pathway")

cg2 <- cg2%>%as.data.frame()%>%
    mutate(Diff.in=as.numeric(gsub("/.*","",GeneRatio)),
           Diff.total=as.numeric(gsub(".*/","",GeneRatio)),
           Diff.not=Diff.total-Diff.in)
cg2 <- cg2%>%
    mutate(Bg.in=as.numeric(gsub("/.*","", BgRatio)),
           Bg.total=as.numeric(gsub(".*/","", BgRatio)),
           Bg.not=Bg.total-Bg.in)

cg2 <- odds.fun(cg2)
cg2 <- cg2%>%dplyr::select(Cluster, p.adjust, odds)%>%full_join(tmp, by=c("Cluster"="rn"))
cg2 <- cg2%>%mutate(p2=ifelse(p.adjust<0.1, 2, 1))
###
fig <- ExampleGOplot(cg2, nbreak=3)+
   ggtitle("Type I IFN signaling")+
   theme(plot.title=element_text(hjust=0.5, size=14))

figfn <-"./4_enrichment.outs/Figure3.2_IFN.png"
png(figfn, width=550, height=400, res=120)
print(fig)
dev.off()


###
###


