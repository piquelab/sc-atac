## TORUS for eQTL
 
We use FastQTL results with 3PCs, which maximize number of egenes, as torus input. 
We also try different kinds of genetic variants annotation. In final, we use the SNP annotation with `pct` 0.02 cell-type active peaks, cell-type specific motifs and response motifs across contrats. The annotation files are in the path `/nfs/rprdata/julong/sc-atac/genetic.analysis_torus_2021-10-14/SNPannnotation/4_SNPAnnot.outs/pct_0.02/3_Bcell_union_torus.annot.gz`. The output are in the path `./torus_output/torus_peak_pct_0.02_union`.
