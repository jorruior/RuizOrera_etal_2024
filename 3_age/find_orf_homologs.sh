#!/bin/bash
#$ -cwd -R yes -l m_mem_free=80G -l h_rt=72:00:00
#$ -o logs/analysis2.o -e logs/analysis2.e

python3 2_find_orf_homologs.py $1 $2 $3
