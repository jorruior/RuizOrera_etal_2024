#!/bin/bash
#$ -cwd -R yes -l m_mem_free=200G -l longrun -l h_rt=200:00:00 -pe smp 1
#$ -o logs/analysis1.o -e logs/analysis1.e

python3 1_find_gene_homologs.py $1 $2 $3 $4
