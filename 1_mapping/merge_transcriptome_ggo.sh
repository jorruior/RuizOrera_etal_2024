#!/bin/bash
#$ -cwd -R yes -l m_mem_free=60G -l h_rt=72:00:00 -pe smp 4
#$ -o logs/merge_bams.o -e logs/merge_bams.e

#source ~/.bash_profile


##Gorilla
stringtie --merge -G annotation/Gorilla_gorilla.gorGor4.98.gtf -m 200 -F 0.5 -f 0.2 -l STRG_GGO -o annotation/Gorilla_gorilla.gorGor4.98.stringtie2021.gtf gg_CM_AssumboCG02___mR_S43_R1_001.fastq.gz/gg_CM_AssumboCG02___mR_S43_R1_001_val_1.fq.gz_stringtie.gtf gg_CM_AssumboDM05___mR_S43_R1_001.fastq.gz/gg_CM_AssumboDM05___mR_S43_R1_001_val_1.fq.gz_stringtie.gtf gg_CM_AssumboSK09___mR_S43_R1_001.fastq.gz/gg_CM_AssumboSK09___mR_S43_R1_001_val_1.fq.gz_stringtie.gtf gg_iPS_Gorilla01__mR_S45_R1_001.fastq.gz/gg_iPS_Gorilla01__mR_S45_R1_001_val_1.fq.gz_stringtie.gtf gg_iPS_Gorilla02__mR_S45_R1_001.fastq.gz/gg_iPS_Gorilla02__mR_S45_R1_001_val_1.fq.gz_stringtie.gtf gg_iPS_Gorilla03__mR_S45_R1_001.fastq.gz/gg_iPS_Gorilla03__mR_S45_R1_001_val_1.fq.gz_stringtie.gtf

gffcompare -R -r annotation/Gorilla_gorilla.gorGor4.98.gtf annotation/Gorilla_gorilla.gorGor4.98.stringtie2021.gtf -o annotation/Gorilla_gorilla.gorGor4.98.stringtie2021.gtf.cuffcomp

python3 ensemblarize_gtf.py annotation/Gorilla_gorilla.gorGor4.98.gtf annotation/Gorilla_gorilla.gorGor4.98.stringtie2021.gtf.cuffcomp.annotated.gtf | sort -k1,1 -k4,4n -k5,5n > annotation/Gorilla_gorilla.gorGor4.98.stringtie2021.fixed.gtf

mkdir annotation/Gorilla_gorilla.gorGor4.dna.toplevel.stringtie2021.fixed

STAR --runThreadN 4 \
--runMode genomeGenerate \
--limitGenomeGenerateRAM=141278166400 \
--genomeDir annotation/Gorilla_gorilla.gorGor4.dna.toplevel.stringtie2021.fixed \
--genomeFastaFiles annotation/Gorilla_gorilla.gorGor4.dna.toplevel.onlychr.fa \
--sjdbGTFfile annotation/Gorilla_gorilla.gorGor4.98.stringtie2021.fixed.gtf \
--sjdbOverhang 75

gffread annotation/Gorilla_gorilla.gorGor4.98.stringtie2021.fixed.gtf -g annotation/Gorilla_gorilla.gorGor4.dna.toplevel.onlychr.fa -w annotation/Gorilla_gorilla.gorGor4.98.stringtie2021.fixed.fa
gffread annotation/Gorilla_gorilla.gorGor4.98.stringtie2021.fixed.gtf -g annotation/Gorilla_gorilla.gorGor4.dna.toplevel.onlychr.fa -x annotation/Gorilla_gorilla.gorGor4.98.stringtie2021.fixed.CDS.fa
python3 translate_orfs.py annotation/Gorilla_gorilla.gorGor4.98.stringtie2021.fixed.CDS.fa > annotation/Gorilla_gorilla.gorGor4.98.stringtie2021.fixed.prot.fa

#Removed chrms: CABD030151492.1, CABD030130064.1

sort -k1,1 -k4,4n -k5,5n annotation/Gorilla_gorilla.gorGor4.98.stringtie2021.fixed.gtf > annotation/Gorilla_gorilla.gorGor4.98.stringtie2021.fixed.sorted.gtf
