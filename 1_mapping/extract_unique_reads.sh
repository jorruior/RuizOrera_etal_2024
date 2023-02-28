#!/bin/bash
#$ -cwd -R yes -l m_mem_free=60G -l h_rt=24:00:00 -pe smp 4
#$ -o logs/mapping2.o -e logs/mapping2.e

#source ~/.bash_profile


cat <(samtools view -H $1) <(samtools view $1 | grep -P "NH:i:1\t") | samtools view -bS - > $1.unique.bam
stringtie $1.unique.bam -M 0.5 -j 3 -p 4 -f 0.2 -G $2 -e -A $1.unique.bam_abundance -o $1.unique.bam_stringtie.gtf
