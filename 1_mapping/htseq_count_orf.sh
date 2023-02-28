#!/bin/bash
#$ -cwd -R yes -l m_mem_free=50G -l h_rt=48:00:00
#$ -o logs/htseq.o -e logs/htseq.e

bam=$1 
gtf=$2
orient=$3 
tag=$4

htseq-count -f bam -s $orient -t CDS --nonunique all -i orf_id $bam $gtf > $bam.$tag.htseq.rna.counts
