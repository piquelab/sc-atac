##
library(tidyverse)
library(data.table)
###
options(scipen=200)

outdir <- "outout"
dir.create(outdir, recursive=T)


### prepare file for liftover-------------------
fn <- "SCAIP1-6_bed.gz"
bed <- fread(fn, header=F, data.table=F)

bed2 <- bed%>%separate(V3, c("chr_pos", "rs"), ";")


######################
### 1.  clean SNP ####
######################

bed_rs <- bed2%>%filter(!is.na(rs), grepl("rs", rs))%>%
    mutate(chr=paste("chr", V1, sep=""), pos=V2, pos2=V2+1)%>%
    dplyr::select(chr, pos, pos2, rs)
##
gfn <- gzfile("SCAIP_clean.bed.gz")
write.table(bed_rs, gfn, sep="\t", row.names=F, col.names=F, quote=F)
## close(gfn)


### run liftOver for clean SNP
### ./liftOver SCAIP_clean.bed.gz hg19ToHg38.over.chain.gz ./output/1.1_clean_output.bed ./output/1.1_clean_unlifted.bed &


####################################################################
### 2. some SNP don't have rs id or id name doesn't include "rs" ### 
####################################################################
bed_na <- bed2%>%filter(!grepl("rs", rs))%>%
 
allele <- str_split(bed_na$chr_pos, ":", simplify=T)

bed_annovar <- data.frame(chr=bed_na$V1, pos=bed_na$V2, pos2=bed_na$V2, ref=allele[,3], alt=allele[,4])
opfn <- "SCAIP_NOrs_annovar.input" 
write.table(bed_annovar, opfn, sep="\t", row.names=F, col.names=F, quote=F)


### run annovar 

### lifover format 
fn <- "SCAIP_NOrs_annovar.input.hg19_avsnp150_dropped"
annovar <- fread(fn, header=F, data.table=F)
##
anvar2 <- annovar%>%filter(grepl("rs", V2))%>%
   mutate(chr=paste("chr", V3, sep=""), pos=V4, pos2=V4+1)%>%
   dplyr::select(chr, pos, pos2, rs=V2) 

### output
gfn <- gzfile("SCAIP_NOrs_recover.bed.gz")
write.table(anvar2, gfn, sep="\t", row.names=F, col.names=F, quote=F)

### run liftover
### ./liftOver SCAIP_NOrs_recover.bed.gz hg19ToHg38.over.chain.gz ./output/2.1_recover_output.bed ./output/2.1_recover_unlifted.bed &


######################
### recover SNP id ###
######################

anvar <- fread("SCAIP_NOrs_annovar.input.hg19_avsnp150_dropped", header=F, data.table=F)
anvar <- anvar%>%mutate(chr_pos=paste(V3, V4, sep="_"))%>%dplyr::select(chr_pos, rs2=V2)

##
bed_new <- bed%>%mutate(chr_pos=paste(V1, V2, sep="_"), rs=bed2$rs)%>%
   dplyr::select(chr=V1, pos=V2, chr_pos, combID=V3, rs)
###
bed_new <- bed_new%>%left_join(anvar, by="chr_pos")%>%
    mutate(rs=ifelse(!grepl("rs",rs), rs2, rs))%>%dplyr::select(-rs2)
###

###
x1 <- fread("./output/1.1_clean_output.bed", header=F, data.table=F)
x2 <- fread("./output/2.1_recover_output.bed", header=F, data.table=F)
x <- rbind(x1, x2)
##
x <- x%>%mutate(chr=gsub("chr", "", V1), chr_pos_grch38=paste(chr, V2, sep="_"))%>%
    dplyr::select(chr_pos_grch38, rs=V4)

bed_new2 <- bed_new%>%drop_na(rs)%>%inner_join(x, by="rs") ## 6,431,934 SNPs

###
gfn <- gzfile("SCAIP_final_bed.gz")
write.table(bed_new2, gfn, row.names=F, col.names=F, quote=F, sep="\t")

bed2 <- fread("SCAIP_final_bed.gz", header=F, data.table=F)

