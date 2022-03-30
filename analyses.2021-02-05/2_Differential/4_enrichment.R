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
theme_set(theme_grey())

rm(list=ls())

###
outdir <- "./4_enrichment.outs/"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)


###
res <- read_rds("./1.2_DiffPeak.outs/2.0_DESeq.results.rds")

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

cg <- read_rds("./4_enrichment.outs/1_enrichGO.rds")

cl <- paste( rep(c("LPS", "PHA", "LPS-DEX", "PHA-DEX"),each=4),
   rep(c("Bcell", "Monocyte", "NKcell", "Tcell"), times=4), sep=".")
cl <- paste(rep(cl,times=2), rep(c("Up","Down"), each=16), sep=".")

ii <- paste(rep(c("A","B","C","D"),each=4),rep(1:4,times=4), sep="")
cl2 <- paste(rep(c("X","Y"),each=16), rep(ii,times=2), sep=".")
cluster2 <- setNames(cl2, cl)
lab2 <- setNames(gsub("-","+",cl),cl2)

##
cg0 <- as.data.frame(cg)
x <- cluster2[as.character(cg0$Cluster)]
cg2 <- cg%>%
   mutate(ClusterNew=x,
          maxGSSize=as.numeric(gsub("/.*", "", BgRatio)),
           ngene=as.numeric(gsub(".*/", "", GeneRatio)))%>%
   dplyr::filter(maxGSSize>10, maxGSSize<500, ngene>10, p.adjust<0.1)

fig1 <- enrichplot::dotplot(cg2, x="ClusterNew", showCategory=5)+
   scale_x_discrete("", labels=lab2)+
   theme(axis.text.x=element_text(angle=60, hjust=1,size=12),
         axis.text.y=element_text(size=10))

### 
figfn <- "./4_enrichment.outs/Figure1.1_GO.png"
png(figfn, width=3500, height=2500, res=160)
print(fig1)
dev.off()


###
figfn <- "./4_enrichment.outs/Figure1.1_GO.pdf"
pdf(figfn, width=20, height=15)
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

cl <- paste( rep(c("LPS", "PHA", "LPS-DEX", "PHA-DEX"),each=4),
   rep(c("Bcell", "Monocyte", "NKcell", "Tcell"), times=4), sep=".")
cl <- paste(rep(cl,times=2), rep(c("Up","Down"), each=16), sep=".")

ii <- paste(rep(c("A","B","C","D"),each=4),rep(1:4,times=4), sep="")
cl2 <- paste(rep(c("X","Y"),each=16), rep(ii,times=2), sep=".")
cluster2 <- setNames(cl2, cl)
lab2 <- setNames(gsub("-","+",cl),cl2)

##
ck0 <- as.data.frame(ck)
x <- cluster2[as.character(ck0$Cluster)]
ck2 <- ck%>%
   mutate(ClusterNew=x,
          maxGSSize=as.numeric(gsub("/.*", "", BgRatio)),
          ngene=as.numeric(gsub(".*/", "", GeneRatio)))%>%
   dplyr::filter(maxGSSize>10, maxGSSize<500, ngene>10, p.adjust<0.1)

fig2 <- enrichplot::dotplot(ck2, x="ClusterNew", showCategory=5)+
   scale_x_discrete("", labels=lab2)+
   theme(axis.text.x=element_text(angle=60, hjust=1,size=12),
         axis.text.y=element_text(size=10))

###
figfn <- "./4_enrichment.outs/Figure2.1_KEGG.png"
png(figfn, width=2000, height=1500, res=150)
print(fig2)
dev.off()

##
figfn <- "./4_enrichment.outs/Figure2.1_KEGG.pdf"
pdf(figfn, width=20, height=15)
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


ExampleGOplot <- function(cg){

### prepare data    
   ## x <- str_split(cg$GeneRatio, "/", simplify=T)
   ## GeneRatio <- as.numeric(x[,1])/as.numeric(x[,2])
   Drt2 <- c("Up"=1, "Down"=2) 
   cg <- cg%>%mutate(Direction2=Drt2[direction.x],
      contrast2=paste(Direction2, contrast.x, sep="."))%>%
      mutate(contrast2=gsub("-", "+", contrast2)) 
   ## cg$size <- rep(1,nrow(cg))
   ## cg$size[GeneRatio>=0.05&GeneRatio<0.15] <- 2
   ## cg$size[GeneRatio>=0.15] <- 3 
   #
   cg <- cg%>%drop_na(odds)
   fig0 <- ggplot(cg, aes(x=contrast2, y=MCls.x))+
      geom_point(aes(size=odds, colour=p2))+
      scale_x_discrete(labels=c("1.LPS"="LPS.Up", "2.LPS"="LPS.Down",
         "1.LPS+DEX"="LPS+DEX.Up", "2.LPS+DEX"="LPS+DEX.Down",
         "1.PHA"="PHA.Up", "2.PHA"="PHA.Down",
         "1.PHA+DEX"="PHA+DEX.Up", "2.PHA+DEX"="PHA+DEX.Down"))+
      scale_colour_gradient(name="p.adjust",                           
         low="blue", high="red", na.value=NA, trans="reverse", n.breaks=5,
         guide=guide_colourbar(order=1))+    #"#ffa500"
      scale_size_binned("odds ratio",
         guide=guide_bins(show.limits=TRUE, axis=TRUE,
           axis.show=arrow(length=unit(1.5,"mm"), ends="both"), order=2),
         n.breaks=4)+
      theme_bw()+
      theme(axis.title=element_blank(),
         ## axis.text.y=element_text(size=12),
         legend.background=element_blank(),
         ## legend.title=element_text(size=8),
         ## legend.text=element_text(size=6),
         legend.key.size=grid::unit(0.6, "lines"))
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
cg2 <- cg2%>%full_join(tmp, by=c("Cluster"="rn"))
cg2 <- cg2%>%mutate(p2=ifelse(p.adjust>0.2, NA, p.adjust))

fig <- ExampleGOplot(cg2)+
   ggtitle("response to lipopolysaccharide")+
   theme(axis.text.x=element_text(angle=-90, size=8, hjust=0, vjust=0.5),
         plot.title=element_text(hjust=0.5, size=12))

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
cg2 <- cg2%>%full_join(tmp, by=c("Cluster"="rn"))
cg2 <- cg2%>%mutate(p2=ifelse(p.adjust>0.2, NA, p.adjust))

fig <- ExampleGOplot(cg2)+
   ggtitle("Type I IFN signaling")+
   theme(axis.text.x=element_text(angle=-90, size=8, hjust=0, vjust=0.5),
         plot.title=element_text(hjust=0.5, size=12))

figfn <-"./4_enrichment.outs/Figure3.2_IFN.png"
png(figfn, width=550, height=400, res=120)
print(fig)
dev.off()


###
###


