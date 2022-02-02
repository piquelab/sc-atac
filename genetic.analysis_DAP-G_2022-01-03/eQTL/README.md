### script for running DAP-G

### 1st step
```ruby
perl batch_process.pl -e ./expressions/Tcell_CTRL.bed.gz -g geno.vcf.gz -c ./covar/Tcell_CTRL_3.GEPCs.txt -t Tcell_CTRL
```

### 2rd step
```ruby
bash Tcell_CTRL.assemble.cmd
```

### 3rd step
