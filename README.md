# sc-atac (Single cell ATAC data analyses)

The working directory is in the path `/nfs/rprdata/julong/sc-atac/`. 
We performed a series of analysis for sc-atac in the following folder. 
- `count.SCAIP.2021-01-14` align reads for each library using `cellranger-atac` from fastq file(`./fastq`)
- `demux.2021-01-23` deconvolute the identity of each cell to determine which individual the cell is from  
- `analyses.2021-02-05` includes rountine analysis of scATAC-seq data, including
  - `1_processing` for preprocessing data, dimensional reduction analysis, cluster analysis, and cell-type annotation 
  - `2_Differential` perform differential chromatin accessibility analysis for pseudo-bulk data 
  - `3_motif` for TF motif occurence matrix, TF motif activity computing and differential motif analysis 
- `genetic.analysis_torus_2021-10-14` 
- `genetic.analysis_DAP-G_2022-01-03`
- `genetic.analysis_ALOFT` contains the scripts for genetic analysis of ALOFT eQTL mapping 
- `genetic.analysis_GTExV8` contains the scripts for genetic analysis of GTEx v8 blood eQTL mapping
- `twas_analysis_2022-10-16` contains the scripts for the identification of asthma risk genes by integration of TWAS and colocalization. 
      
   