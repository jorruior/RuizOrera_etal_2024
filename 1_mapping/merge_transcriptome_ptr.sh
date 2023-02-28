#!/bin/bash
#$ -cwd -R yes -l m_mem_free=60G -l h_rt=72:00:00 -pe smp 4
#$ -o logs/merge_bams.o -e logs/merge_bams.e

#source ~/.bash_profile


##Chimp
stringtie --merge -G annotation/Pan_troglodytes.Pan_tro_3.0.98.gtf -m 200 -F 0.5 -f 0.2 -l STRG_PTR -o annotation/Pan_troglodytes.Pan_tro_3.0.98.stringtie2021.gtf pt_CM_KarlchenCG23__mR_S44_R1_001.fastq.gz/pt_CM_KarlchenCG23__mR_S44_R1_001_val_1.fq.gz_stringtie.gtf pt_CM_KarlchenDM02__mR_S44_R1_001.fastq.gz/pt_CM_KarlchenDM02__mR_S44_R1_001_val_1.fq.gz_stringtie.gtf pt_CM_KarlchenDM07__mR_S44_R1_001.fastq.gz/pt_CM_KarlchenDM07__mR_S44_R1_001_val_1.fq.gz_stringtie.gtf pt_lv_Barney____mR_S37_R1_001.fastq.gz/pt_lv_Barney____mR_S37_R1_001_val_1.fq.gz_stringtie.gtf pt_lv_Dirk1_____mR_S35_R1_001.fastq.gz/pt_lv_Dirk1_____mR_S35_R1_001_val_1.fq.gz_stringtie.gtf pt_lv_Emmanuel__mR_S34_R1_001.fastq.gz/pt_lv_Emmanuel__mR_S34_R1_001_val_1.fq.gz_stringtie.gtf pt_lv_Fangina____mR_S36_R1_001.fastq.gz/pt_lv_Fangina____mR_S36_R1_001_val_1.fq.gz_stringtie.gtf pt_lv_Theo______mR_S38_R1_001.fastq.gz/pt_lv_Theo______mR_S38_R1_001_val_1.fq.gz_stringtie.gtf pt_iPS_Chimp01____mR_S46_R1_001.fastq.gz/pt_iPS_Chimp01____mR_S46_R1_001_val_1.fq.gz_stringtie.gtf pt_iPS_Chimp02____mR_S46_R1_001.fastq.gz/pt_iPS_Chimp02____mR_S46_R1_001_val_1.fq.gz_stringtie.gtf pt_iPS_Chimp03____mR_S46_R1_001.fastq.gz/pt_iPS_Chimp03____mR_S46_R1_001_val_1.fq.gz_stringtie.gtf

gffcompare -R -r annotation/Pan_troglodytes.Pan_tro_3.0.98.gtf annotation/Pan_troglodytes.Pan_tro_3.0.98.stringtie2021.gtf -o annotation/Pan_troglodytes.Pan_tro_3.0.98.stringtie2021.gtf.cuffcomp

python3 ensemblarize_gtf.py annotation/Pan_troglodytes.Pan_tro_3.0.98.gtf annotation/Pan_troglodytes.Pan_tro_3.0.98.stringtie2021.gtf.cuffcomp.annotated.gtf | sort -k1,1 -k4,4n -k5,5n > annotation/Pan_troglodytes.Pan_tro_3.0.98.stringtie2021.fixed.gtf

mkdir annotation/Pan_troglodytes.Pan_tro_3.0.dna.toplevel.stringtie2021.fixed

STAR --runThreadN 4 \
--runMode genomeGenerate \
--limitGenomeGenerateRAM=141278166400 \
--genomeDir annotation/Pan_troglodytes.Pan_tro_3.0.dna.toplevel.stringtie2021.fixed \
--genomeFastaFiles annotation/Pan_troglodytes.Pan_tro_3.0.dna.toplevel.onlychr.fa \
--sjdbGTFfile annotation/Pan_troglodytes.Pan_tro_3.0.98.stringtie2021.fixed.gtf \
--sjdbOverhang 75

gffread annotation/Pan_troglodytes.Pan_tro_3.0.98.stringtie2021.fixed.gtf -g annotation/Pan_troglodytes.Pan_tro_3.0.dna.toplevel.onlychr.fa -w annotation/Pan_troglodytes.Pan_tro_3.0.98.stringtie2021.fixed.fa
gffread annotation/Pan_troglodytes.Pan_tro_3.0.98.stringtie2021.fixed.gtf -g annotation/Pan_troglodytes.Pan_tro_3.0.dna.toplevel.onlychr.fa -x annotation/Pan_troglodytes.Pan_tro_3.0.98.stringtie2021.fixed.CDS.fa
python3 translate_orfs.py annotation/Pan_troglodytes.Pan_tro_3.0.98.stringtie2021.fixed.CDS.fa > annotation/Pan_troglodytes.Pan_tro_3.0.98.stringtie2021.fixed.prot.fa

sort -k1,1 -k4,4n -k5,5n annotation/Pan_troglodytes.Pan_tro_3.0.98.stringtie2021.fixed.gtf > annotation/Pan_troglodytes.Pan_tro_3.0.98.stringtie2021.fixed.sorted.gtf
