#!/bin/bash
#$ -cwd -R yes -l m_mem_free=60G -l h_rt=72:00:00 -pe smp 4
#$ -o logs/merge_bams.o -e logs/merge_bams.e

#source ~/.bash_profile


##Macaque
stringtie --merge -G annotation/Macaca_mulatta.Mmul_10.98.gtf -m 200 -F 0.5 -f 0.2 -l STRG_MML -o annotation/Macaca_mulatta.Mmul_10.98.stringtie2021.gtf m_CM_WT101______mR_S01_R1_001.fastq.gz/rm_CM_WT101______mR_S01_R1_001_val_1.fq.gz_stringtie.gtf rm_CM_WT102______mR_S01_R1_001.fastq.gz/rm_CM_WT102______mR_S01_R1_001_val_1.fq.gz_stringtie.gtf rm_CM_WT103______mR_S01_R1_001.fastq.gz/rm_CM_WT103______mR_S01_R1_001_val_1.fq.gz_stringtie.gtf rm_lv_Ri9809011_mR_S39_R1_001.fastq.gz/rm_lv_Ri9809011_mR_S39_R1_001_val_1.fq.gz_stringtie.gtf rm_lv_RiN8______mR_S41_R1_001.fastq.gz/rm_lv_RiN8______mR_S41_R1_001_val_1.fq.gz_stringtie.gtf rm_iPS_WT101_____mR_S01_R1_001.fastq.gz/rm_iPS_WT101_____mR_S01_R1_001_val_1.fq.gz_stringtie.gtf rm_iPS_WT102_____mR_S01_R1_001.fastq.gz/rm_iPS_WT102_____mR_S01_R1_001_val_1.fq.gz_stringtie.gtf rm_iPS_WT103_____mR_S01_R1_001.fastq.gz/rm_iPS_WT103_____mR_S01_R1_001_val_1.fq.gz_stringtie.gtf

gffcompare -R -r annotation/Macaca_mulatta.Mmul_10.98.gtf annotation/Macaca_mulatta.Mmul_10.98.stringtie2021.gtf -o annotation/Macaca_mulatta.Mmul_10.98.stringtie2021.gtf.cuffcomp

python3 ensemblarize_gtf.py annotation/Macaca_mulatta.Mmul_10.98.gtf annotation/Macaca_mulatta.Mmul_10.98.stringtie2021.gtf.cuffcomp.annotated.gtf | sort -k1,1 -k4,4n -k5,5n > annotation/Macaca_mulatta.Mmul_10.98.stringtie2021.fixed.gtf

mkdir annotation/Macaca_mulatta.Mmul_10.dna.toplevel.stringtie2021.fixed

STAR --runThreadN 4 \
--runMode genomeGenerate \
--limitGenomeGenerateRAM=141278166400 \
--genomeDir annotation/Macaca_mulatta.Mmul_10.dna.toplevel.stringtie2021.fixed \
--genomeFastaFiles annotation/Macaca_mulatta.Mmul_10.dna.toplevel.fa \
--sjdbGTFfile annotation/Macaca_mulatta.Mmul_10.98.stringtie2021.fixed.gtf \
--sjdbOverhang 75

gffread annotation/Macaca_mulatta.Mmul_10.98.stringtie2021.fixed.gtf -g annotation/Macaca_mulatta.Mmul_10.dna.toplevel.fa -w annotation/Macaca_mulatta.Mmul_10.98.stringtie2021.fixed.fa
gffread annotation/Macaca_mulatta.Mmul_10.98.stringtie2021.fixed.gtf -g annotation/Macaca_mulatta.Mmul_10.dna.toplevel.fa -x annotation/Macaca_mulatta.Mmul_10.98.stringtie2021.fixed.CDS.fa

sort -k1,1 -k4,4n -k5,5n annotation/Macaca_mulatta.Mmul_10.98.stringtie2021.fixed.gtf > annotation/Macaca_mulatta.Mmul_10.98.stringtie2021.fixed.sorted.gtf
