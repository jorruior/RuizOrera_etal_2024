#!/bin/bash
#$ -cwd -R yes -l m_mem_free=120G -l h_rt=50:00:00
#$ -o logs/orfmapper.o -e logs/orfmapper.e

python3 consensus_ORFs.py ../1_mapping/$1\_pooled/ribo_pooled.bam.collapsed lists/$1.txt X X X X 0 0 $1 only_collapse

python3 ORF_mapper_to_GENCODE_primates.py -C yes -m psite_overlap -s ../1_mapping/$1\_pooled/rna_pooled.stringtie.gtf -d $2.stringtie2021 -f ../1_mapping/$1\_pooled/ribo_pooled.collapsed.fa.collapsed -b ../1_mapping/$1\_pooled/ribo_pooled.collapsed.bed.collapsed -o ../1_mapping/$1\_pooled/$1 -l 1

python3 consensus_ORFs.py ../1_mapping/$1\_pooled/$1.orfs.gtf lists/$1.txt ../1_mapping/$1\_pooled/$1.orfs.fa ../1_mapping/$1\_pooled/$1.orfs.allframes.bed none ../1_mapping/annotation/$2.gtf 2 0 $1 no

