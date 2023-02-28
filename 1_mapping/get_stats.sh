#!/bin/bash

#1. Folder with files
#2. Output name

source ~/.bash_profile

cd $1
mkdir $2
#multiqc .
cp *html $2/
cp */*html $2/

#Filtering stats
for f in `ls */*trimmed.fq.gz | grep -v val_1 | grep _sec`; do zcat $f | wc -l | awk '{print $1/4}'; done > $2/trimmed_stats.txt

#Mapped stats 
for f in `ls */*gz_STAR*/Log.final.out | grep -v val_1 | grep _sec`; do sed -n -e 6p -e 9p -e 10p $f | tr '\n' '\t' | awk -v a=$f 'BEGIN{FS="\t"}{print a"\t"$2"\t"$4"\t"$6}' | sed 's#Log.final.out##g'; done > $2/map_stats.txt
echo -e "sample\ttrimmed_reads\tfiltered_reads\tmapped_reads\tmapped_%" > $2/total_stats.txt
paste $2/map_stats.txt $2/trimmed_stats.txt | awk '{print $1"\t"$5"\t"$2"\t"$3"\t"$4}' >> $2/total_stats.txt
#rm $2/map_stats.txt
#rm $2/trimmed_stats.txt

#Periodicity stats
echo -e "sample\tperiodicities\tread_lengths\t%_reads" >> $2/periodicity_stats.txt
for f in */*calcs; do sort -nrk10,10 $f | grep -v FALSE | grep TRUE | awk '{print int($2+0.5)"\t"$3"\t"$10}' | python -c "import sys; print('\n'.join(','.join(c) for c in zip(*(l.split() for l in sys.stdin.readlines() if l.strip()))))" | tr '\n' '\t' | awk -v a=$f '{print a"\t"$0}'; done >> $2/periodicity_stats.txt

