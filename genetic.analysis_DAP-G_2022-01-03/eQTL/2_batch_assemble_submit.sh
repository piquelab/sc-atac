#/bin/bash/

cat datasets.txt | \
while read condition;
do
   echo ${condition}
   sbatch -q express -N 1-1 -n 1 --mem=10G --time=12:00:00 --job-name=${condition} --wrap "bash ${condition}.assemble.cmd"
   sleep 1;
   # chmod a+x ${condition}.assemble.cmd
done

