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
          peak_region=paste(chr,start,end,sep="-"))%>%
   dplyr::select(peak_region, geneId, distanceToTSS)# flank_geneIds) 
    

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
    mutate(comb=paste(MCls, contrast, sep="_"))%>%
    dplyr::filter(abs(estimate)>0.5, p.adjusted<0.1, MCls!="DC")
res2 <- res2%>%left_join(peakAnno, by="peak_region")

##
comb <- sort(unique(res2$comb))
df0 <- map_dfr(comb, function(ii){
   tmp <- res2%>%dplyr::filter(comb==ii)%>%mutate(dis=abs(distanceToTSS))%>%
      group_by(geneId)%>%
      dplyr::filter(dis==min(dis,na.rm=T))%>%ungroup()
   tmp
})

df0 <- as.data.frame(df0)

x <- bitr(df0$geneId, fromType="ENSEMBL",
          toType=c("ENTREZID","SYMBOL"), OrgDb=org.Hs.eg.db)

geneCluster <- df0%>%inner_join(x, by=c("geneId"="ENSEMBL"))%>%
   mutate(direction=ifelse(estimate>0, "Up", "Down"))


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
          maxGSSize=as.numeric(gsub("/.*", "", BgRatio)))%>%
   dplyr::filter(maxGSSize<500, p.adjust<0.1)

fig1 <- enrichplot::dotplot(cg2, x="ClusterNew", showCategory=5)+
   scale_x_discrete("", labels=lab2)+
   theme(axis.text.x=element_text(angle=60, hjust=1,size=12),
         axis.text.y=element_text(size=10))

###
figfn <- "./4_enrichment.outs/Figure1.1_GO.png"
png(figfn, width=2000, height=1600, res=120)
print(fig1)
dev.off()


###
figfn <- "./4_enrichment.outs/Figure1.2_GO.pdf"
pdf(figfn, width=20, height=15)
print(fig1)
dev.off()















