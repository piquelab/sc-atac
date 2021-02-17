library(tidyverse)
library(purrr)
library(ggplot2)
library(cowplot)
theme_set(theme_grey())
##bash script for summary results 
## grep 'anno' ./*/outs/summary.csv | head -1 | sed 's/.*csv:/Library ID,/' > all_summary.csv;
## grep -v 'anno' ./*/outs/summary.csv | sed 's/\/outs.*csv:/,/;s/\.\///' >> all_summary.csv
##ls ./*/outs/summary.csv
output <- "./output/"
if (!file.exists(output)) dir.create(output,showWarnings=F)

###
folders <- dir(".", "^SCAIP6.*")
aa <- map_dfr(folders, function(x){
  fn<- paste(x, "/outs/summary.csv", sep="")
  a <- read.csv(fn)%>%mutate(EXP=x)
})


###figure 1
fig1 <- ggplot(aa,aes(x=EXP, y=num_fragments))+
        geom_bar(stat="identity", fill="#1c9099")+
        ggtitle("total number of read pairs")+
        theme_bw()+
        theme(axis.text.x=element_text(angle=60, hjust=1, size=8),
              axis.title=element_blank(),
              plot.title=element_text(hjust=0.5))

figfn <- "./output/Figure1.total_pairs.png"
png(figfn, width=500, height=500, res=120)
print(fig1)
dev.off() 

### numbers of cells
fig2 <- ggplot(aa,aes(x=EXP, y=annotated_cells))+
        geom_bar(stat="identity", fill="#1c9099")+
        ggtitle("Number of cells")+
        theme_bw()+
        theme(axis.text.x=element_text(angle=60, hjust=1, size=8),
              axis.title=element_blank(),
              plot.title=element_text(hjust=0.5))

figfn <- "./output/Figure2.cells.png"
png(figfn, width=500, height=500, res=120)
print(fig2)
dev.off()

###
fig3 <- ggplot(aa,aes(x=EXP, y=median_fragments_per_cell))+
        geom_bar(stat="identity", fill="#1c9099")+
        ggtitle("Median fragments per cell")+
        theme_bw()+
        theme(axis.text.x=element_text(angle=60, hjust=1, size=8),
              axis.title=element_blank(),
              plot.title=element_text(hjust=0.5))

figfn <- "./output/Figure3.median.fragments.png"
png(figfn, width=500, height=500, res=120)
print(fig3)
dev.off()    
   
###4
fig4 <- ggplot(aa,aes(x=EXP, y=tss_enrichment_score))+
        geom_bar(stat="identity", fill="#1c9099")+
        ggtitle("tss enrichment score")+
        theme_bw()+
        theme(axis.text.x=element_text(angle=60, hjust=1, size=8),
              axis.title=element_blank(),
              plot.title=element_text(hjust=0.5))

figfn <- "./output/Figure4.tss_score.png"
png(figfn, width=500, height=500, res=120)
print(fig4)
dev.off() 



         