### Documents for SNP annotation

### 1st step

### 2rd step

### 3rd step
```ruby
cd ./Motif_file/
split motifList2020.txt -l 30 -d -a 3 splitMotif
ls | grep `^split` > ../motif_split.txt
```

### annotation files
Default using peaks expressed at least 2% cells 
- cell-type specific+union of response motifs `./4_SNPAnnot.outs/pct_0.02/3_*_union_torus.annot.gz`
- cell-type specific peaks annotation `./4_SNPAnnot.outs/pct_0.02/*_Active_peak_torus.annot.gz`

### annotation motif each by each motif
The path is in `./5_SNPAnnot.outs/pct_0.02_*/`

We have three types of annotation files for different purposes, 
- `./4_SNPAnnot.outs/pct_0.02/3_*_union_torus.annot.gz` these are cell type specific annotation files with three simplified category,  cell type active peaks, cell type active motifs and response motifs in specific cell type.   
- `./4_SNPAnnot.outs/pct_0.02_combineNew/combine2_torus.annot.gz` is the annotation file with multiple columns, one for any cell type active peak, four columns for cell type active motif and 16 columns for 16 response motifs.
- `./4_SNPAnnot.outs/pct_0.02_combineNew/Union_torus.annot.gz` is the simplified annotation with 1 column, inlcuding any cell type peak, cell type active motif and response motif.

    
  

