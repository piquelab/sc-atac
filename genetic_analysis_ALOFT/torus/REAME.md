# torus for ALOFT


### ALOFT_eQTL results `/nfs/rprdata/ALOFT/AL1-6_ln/FastQTL_corrected/nominals/output/PC1-18.nominals.eQTL.txt.gz`
### FastQTL script `/nfs/rprdata/ALOFT/AL1-6_ln/FastQTL_corrected/nominals/run.FastQTL.txt` default 1MB 

The file structure of directory
- The `1_prepare_map.R` file used for preparing torus input files;
- The `2.5_submit.sh` file used for submitting jobs for running torus program;
- The `2.5_summary.R` file used for summarizing output of torus.  

The annotation files in the directory `./torus_annot/`
- The `combine_torus.annot.gz` file includes 9 annotation, 1 peaking_d, 4 cell-type motif and 4 response motif;
- The `combine2_torus.annot.gz` file includes 17 annotation, 1 peaking_d, 16 condition motif annotation;
- The `Union_torus.annot.gz` file includes 1 annotation.  
