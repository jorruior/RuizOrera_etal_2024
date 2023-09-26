#!/bin/bash
#$ -cwd -R yes -l m_mem_free=350G -l h_rt=48:00:00
#$ -o logs/analysis3.o -e logs/analysis3.e

python3 3_find_region_homologs.py $1 $2 $3 $4
