#!/bin/bash
#$ -cwd -R yes -l m_mem_free=60G -l h_rt=72:00:00 -pe smp 4
#$ -o logs/merge_bams.o -e logs/merge_bams.e

#source ~/.bash_profile


##Human
stringtie --merge -G annotation/Homo_sapiens.GRCh38.98.gtf -m 200 -F 0.5 -f 0.2 -l STRG_HSA -o annotation/Homo_sapiens.GRCh38.98.stringtie2021.gtf hs_CM_BIHi242-C1_mR_S01_R1_001.fastq.gz/hs_CM_BIHi242-C1_mR_S01_R1_001_val_1.fq.gz_stringtie.gtf hs_CM_BIHi242-C2_mR_S01_R1_001.fastq.gz/hs_CM_BIHi242-C2_mR_S01_R1_001_val_1.fq.gz_stringtie.gtf hs_CM_BIHi242-C3_mR_S01_R1_001.fastq.gz/hs_CM_BIHi242-C3_mR_S01_R1_001_val_1.fq.gz_stringtie.gtf hs_CM_BIHi242-C4_mR_S01_R1_001.fastq.gz/hs_CM_BIHi242-C4_mR_S01_R1_001_val_1.fq.gz_stringtie.gtf hs_CM_BIHi242-C5_mR_S01_R1_001.fastq.gz/hs_CM_BIHi242-C5_mR_S01_R1_001_val_1.fq.gz_stringtie.gtf hs_lv_066_C_mR_R1_raw.fastq.gz/hs_lv_066_C_mR_R1_raw_val_1.fq.gz_stringtie.gtf hs_lv_067_C_mR_R1_raw.fastq.gz/hs_lv_067_C_mR_R1_raw_val_1.fq.gz_stringtie.gtf hs_lv_068_C_mR_R1_raw.fastq.gz/hs_lv_068_C_mR_R1_raw_val_1.fq.gz_stringtie.gtf hs_lv_069_C_mR_R1_raw.fastq.gz/hs_lv_069_C_mR_R1_raw_val_1.fq.gz_stringtie.gtf hs_lv_070_C_mR_R1_raw.fastq.gz/hs_lv_070_C_mR_R1_raw_val_1.fq.gz_stringtie.gtf hs_lv_071_C_mR_R1_raw.fastq.gz/hs_lv_071_C_mR_R1_raw_val_1.fq.gz_stringtie.gtf hs_lv_072_C_mR_R1_raw.fastq.gz/hs_lv_072_C_mR_R1_raw_val_1.fq.gz_stringtie.gtf hs_lv_073_C_mR_R1_raw.fastq.gz/hs_lv_073_C_mR_R1_raw_val_1.fq.gz_stringtie.gtf hs_lv_074_C_mR_R1_raw.fastq.gz/hs_lv_074_C_mR_R1_raw_val_1.fq.gz_stringtie.gtf hs_lv_075_C_mR_R1_raw.fastq.gz/hs_lv_075_C_mR_R1_raw_val_1.fq.gz_stringtie.gtf hs_lv_076_C_mR_R1_raw.fastq.gz/hs_lv_076_C_mR_R1_raw_val_1.fq.gz_stringtie.gtf hs_lv_077_C_mR_R1_raw.fastq.gz/hs_lv_077_C_mR_R1_raw_val_1.fq.gz_stringtie.gtf hs_lv_078_C_mR_R1_raw.fastq.gz/hs_lv_078_C_mR_R1_raw_val_1.fq.gz_stringtie.gtf hs_lv_079_C_mR_R1_raw.fastq.gz/hs_lv_079_C_mR_R1_raw_val_1.fq.gz_stringtie.gtf hs_lv_081_C_mR_R1_raw.fastq.gz/hs_lv_081_C_mR_R1_raw_val_1.fq.gz_stringtie.gtf hs_iPS_BIHi242-C01_mR_S01_R1_001.fastq.gz/hs_iPS_BIHi242-C01_mR_S01_R1_001_val_1.fq.gz_stringtie.gtf hs_iPS_BIHi242-C02_mR_S01_R1_001.fastq.gz/hs_iPS_BIHi242-C02_mR_S01_R1_001_val_1.fq.gz_stringtie.gtf hs_iPS_BIHi242-C03_mR_S01_R1_001.fastq.gz/hs_iPS_BIHi242-C03_mR_S01_R1_001_val_1.fq.gz_stringtie.gtf

gffcompare -R -r annotation/Homo_sapiens.GRCh38.98.gtf annotation/Homo_sapiens.GRCh38.98.stringtie2021.gtf -o annotation/Homo_sapiens.GRCh38.98.stringtie2021.gtf.cuffcomp

python3 ensemblarize_gtf.py annotation/Homo_sapiens.GRCh38.98.gtf annotation/Homo_sapiens.GRCh38.98.stringtie2021.gtf.cuffcomp.annotated.gtf | sort -k1,1 -k4,4n -k5,5n > annotation/Homo_sapiens.GRCh38.98.stringtie2021.fixed.gtf

mkdir annotation/Homo_sapiens.GRCh38.dna.toplevel.stringtie2021.fixed 

STAR --runThreadN 4 \
--runMode genomeGenerate \
--limitGenomeGenerateRAM=141278166400 \
--genomeDir annotation/Homo_sapiens.GRCh38.dna.toplevel.stringtie2021.fixed \
--genomeFastaFiles annotation/Homo_sapiens.GRCh38.dna.toplevel.fa \
--sjdbGTFfile annotation/Homo_sapiens.GRCh38.98.stringtie2021.fixed.gtf \
--sjdbOverhang 75

gffread annotation/Homo_sapiens.GRCh38.98.stringtie2021.fixed.gtf -g annotation/Homo_sapiens.GRCh38.dna.toplevel.fa -w annotation/Homo_sapiens.GRCh38.98.stringtie2021.fixed.fa
gffread annotation/Homo_sapiens.GRCh38.98.stringtie2021.fixed.gtf -g annotation/Homo_sapiens.GRCh38.dna.toplevel.fa -x annotation/Homo_sapiens.GRCh38.98.stringtie2021.fixed.CDS.fa
python3 translate_orfs.py annotation/Homo_sapiens.GRCh38.98.stringtie2021.fixed.CDS.fa > annotation/Homo_sapiens.GRCh38.98.stringtie2021.fixed.prot.fa

sort -k1,1 -k4,4n -k5,5n annotation/Homo_sapiens.GRCh38.98.stringtie2021.fixed.gtf > annotation/Homo_sapiens.GRCh38.98.stringtie2021.fixed.sorted.gtf
