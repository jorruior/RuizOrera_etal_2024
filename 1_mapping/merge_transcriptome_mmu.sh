#!/bin/bash
#$ -cwd -R yes -l m_mem_free=60G -l h_rt=72:00:00 -pe smp 4
#$ -o logs/merge_bams.o -e logs/merge_bams.e

#source ~/.bash_profile


##Human
stringtie --merge -G annotation/Mus_musculus.GRCm38.98.gtf -m 200 -F 0.5 -f 0.2 -l STRG_MMU -o annotation/Mus_musculus.GRCm38.98.stringtie2021.gtf 1_mapping/mm_heart_001_mR_R1_raw.fastq.gz/mm_heart_001_mR_R1_raw_val_1.fq.gz_stringtie.gtf 1_mapping/mm_heart_002_mR_R1_raw.fastq.gz/mm_heart_002_mR_R1_raw_val_1.fq.gz_stringtie.gtf 1_mapping/mm_heart_003_mR_R1_raw.fastq.gz/mm_heart_003_mR_R1_raw_val_1.fq.gz_stringtie.gtf 1_mapping/mm_heart_004_mR_R1_raw.fastq.gz/mm_heart_004_mR_R1_raw_val_1.fq.gz_stringtie.gtf 1_mapping/mm_heart_006_mR_R1_raw.fastq.gz/mm_heart_006_mR_R1_raw_val_1.fq.gz_stringtie.gtf 1_mapping/mm_heart_005_mR_R1_raw.fastq.gz/mm_heart_005_mR_R1_raw_val_1.fq.gz_stringtie.gtf 1742sTS.Mouse.Heart.0dpb.Female.fastq.gz/1742sTS.Mouse.Heart.0dpb.Female_trimmed.fq.gz_stringtie.gtf 1902sTS.Mouse.Heart.0dpb.Female.fastq.gz/1902sTS.Mouse.Heart.0dpb.Female_trimmed.fq.gz_stringtie.gtf 1908sTS.Mouse.Heart.0dpb.Male.fastq.gz/1908sTS.Mouse.Heart.0dpb.Male_trimmed.fq.gz_stringtie.gtf 1914sTS.Mouse.Heart.0dpb.Male.fastq.gz/1914sTS.Mouse.Heart.0dpb.Male_trimmed.fq.gz_stringtie.gtf 1920sTS.Mouse.Heart.3dpb.Female.fastq.gz/1920sTS.Mouse.Heart.3dpb.Female_trimmed.fq.gz_stringtie.gtf 1926sTS.Mouse.Heart.3dpb.Male.fastq.gz/1926sTS.Mouse.Heart.3dpb.Male_trimmed.fq.gz_stringtie.gtf 1932sTS.Mouse.Heart.3dpb.Female.fastq.gz/1932sTS.Mouse.Heart.3dpb.Female_trimmed.fq.gz_stringtie.gtf 1950sTS.Mouse.Heart.9wpb.Female.fastq.gz/1950sTS.Mouse.Heart.9wpb.Female_trimmed.fq.gz_stringtie.gtf 1956sTS.Mouse.Heart.9wpb.Male.fastq.gz/1956sTS.Mouse.Heart.9wpb.Male_trimmed.fq.gz_stringtie.gtf 2459sTS.Mouse.Heart.14.5.Male.fastq.gz/2459sTS.Mouse.Heart.14.5.Male_trimmed.fq.gz_stringtie.gtf 2461sTS.Mouse.Heart.14.5.Male.fastq.gz/2461sTS.Mouse.Heart.14.5.Male_trimmed.fq.gz_stringtie.gtf 2487sTS.Mouse.Heart.14.5.Female.fastq.gz/2487sTS.Mouse.Heart.14.5.Female_trimmed.fq.gz_stringtie.gtf 2495sTS.Mouse.Heart.10.5.Male.fastq.gz/2495sTS.Mouse.Heart.10.5.Male_trimmed.fq.gz_stringtie.gtf 2504sTS.Mouse.Heart.15.5.Male.fastq.gz/2504sTS.Mouse.Heart.15.5.Male_trimmed.fq.gz_stringtie.gtf 2506sTS.Mouse.Heart.15.5.Male.fastq.gz/2506sTS.Mouse.Heart.15.5.Male_trimmed.fq.gz_stringtie.gtf 2508sTS.Mouse.Heart.15.5.Female.fastq.gz/2508sTS.Mouse.Heart.15.5.Female_trimmed.fq.gz_stringtie.gtf 2510sTS.Mouse.Heart.15.5.Female.fastq.gz/2510sTS.Mouse.Heart.15.5.Female_trimmed.fq.gz_stringtie.gtf 2512sTS.Mouse.Heart.16.5.Male.fastq.gz/2512sTS.Mouse.Heart.16.5.Male_trimmed.fq.gz_stringtie.gtf 2514sTS.Mouse.Heart.16.5.Male.fastq.gz/2514sTS.Mouse.Heart.16.5.Male_trimmed.fq.gz_stringtie.gtf 2581sTS.Mouse.Heart.10.5.Female.fastq.gz/2581sTS.Mouse.Heart.10.5.Female_trimmed.fq.gz_stringtie.gtf 2585sTS.Mouse.Heart.16.5.Female.fastq.gz/2585sTS.Mouse.Heart.16.5.Female_trimmed.fq.gz_stringtie.gtf 2587sTS.Mouse.Heart.16.5.Female.fastq.gz/2587sTS.Mouse.Heart.16.5.Female_trimmed.fq.gz_stringtie.gtf 2589sTS.Mouse.Heart.17.5.Male.fastq.gz/2589sTS.Mouse.Heart.17.5.Male_trimmed.fq.gz_stringtie.gtf 2591sTS.Mouse.Heart.17.5.Male.fastq.gz/2591sTS.Mouse.Heart.17.5.Male_trimmed.fq.gz_stringtie.gtf 2593sTS.Mouse.Heart.17.5.Female.fastq.gz/2593sTS.Mouse.Heart.17.5.Female_trimmed.fq.gz_stringtie.gtf 2595sTS.Mouse.Heart.17.5.Female.fastq.gz/2595sTS.Mouse.Heart.17.5.Female_trimmed.fq.gz_stringtie.gtf 2622sTS.Mouse.Heart.11.5.Male.fastq.gz/2622sTS.Mouse.Heart.11.5.Male_trimmed.fq.gz_stringtie.gtf 2633sTS.Mouse.Heart.18.5.Male.fastq.gz/2633sTS.Mouse.Heart.18.5.Male_trimmed.fq.gz_stringtie.gtf 2635sTS.Mouse.Heart.18.5.Male.fastq.gz/2635sTS.Mouse.Heart.18.5.Male_trimmed.fq.gz_stringtie.gtf 2639sTS.Mouse.Heart.18.5.Female.fastq.gz/2639sTS.Mouse.Heart.18.5.Female_trimmed.fq.gz_stringtie.gtf 2641sTS.Mouse.Heart.3dpb.Male.fastq.gz/2641sTS.Mouse.Heart.3dpb.Male_trimmed.fq.gz_stringtie.gtf 2643sTS.Mouse.Heart.2wpb.Male.fastq.gz/2643sTS.Mouse.Heart.2wpb.Male_trimmed.fq.gz_stringtie.gtf 2645sTS.Mouse.Heart.2wpb.Male.fastq.gz/2645sTS.Mouse.Heart.2wpb.Male_trimmed.fq.gz_stringtie.gtf 2647sTS.Mouse.Heart.10.5.Female.fastq.gz/2647sTS.Mouse.Heart.10.5.Female_trimmed.fq.gz_stringtie.gtf 2651sTS.Mouse.Heart.2wpb.Female.fastq.gz/2651sTS.Mouse.Heart.2wpb.Female_trimmed.fq.gz_stringtie.gtf 2653sTS.Mouse.Heart.2wpb.Female.fastq.gz/2653sTS.Mouse.Heart.2wpb.Female_trimmed.fq.gz_stringtie.gtf 2655sTS.Mouse.Heart.4wpb.Female.fastq.gz/2655sTS.Mouse.Heart.4wpb.Female_trimmed.fq.gz_stringtie.gtf 2657sTS.Mouse.Heart.4wpb.Female.fastq.gz/2657sTS.Mouse.Heart.4wpb.Female_trimmed.fq.gz_stringtie.gtf 2659sTS.Mouse.Heart.4wpb.Male.fastq.gz/2659sTS.Mouse.Heart.4wpb.Male_trimmed.fq.gz_stringtie.gtf 2661sTS.Mouse.Heart.4wpb.Male.fastq.gz/2661sTS.Mouse.Heart.4wpb.Male_trimmed.fq.gz_stringtie.gtf 2663sTS.Mouse.Heart.9wpb.Male.fastq.gz/2663sTS.Mouse.Heart.9wpb.Male_trimmed.fq.gz_stringtie.gtf 2669sTS.Mouse.Heart.9wpb.Female.fastq.gz/2669sTS.Mouse.Heart.9wpb.Female_trimmed.fq.gz_stringtie.gtf 2815sTS.Mouse.Heart.13.5.Female.fastq.gz/2815sTS.Mouse.Heart.13.5.Female_trimmed.fq.gz_stringtie.gtf 2823sTS.Mouse.Heart.13.5.Male.fastq.gz/2823sTS.Mouse.Heart.13.5.Male_trimmed.fq.gz_stringtie.gtf 2831sTS.Mouse.Heart.13.5.Female.fastq.gz/2831sTS.Mouse.Heart.13.5.Female_trimmed.fq.gz_stringtie.gtf 2839sTS.Mouse.Heart.13.5.Male.fastq.gz/2839sTS.Mouse.Heart.13.5.Male_trimmed.fq.gz_stringtie.gtf 2938sTS.Mouse.Heart.11.5.Female.fastq.gz/2938sTS.Mouse.Heart.11.5.Female_trimmed.fq.gz_stringtie.gtf 2944sTS.Mouse.Heart.11.5.Male.fastq.gz/2944sTS.Mouse.Heart.11.5.Male_trimmed.fq.gz_stringtie.gtf 3031sTS.Mouse.Heart.11.5.Female.fastq.gz/3031sTS.Mouse.Heart.11.5.Female_trimmed.fq.gz_stringtie.gtf 3036sTS.Mouse.Heart.12.5.Male.fastq.gz/3036sTS.Mouse.Heart.12.5.Male_trimmed.fq.gz_stringtie.gtf 3041sTS.Mouse.Heart.12.5.Male.fastq.gz/3041sTS.Mouse.Heart.12.5.Male_trimmed.fq.gz_stringtie.gtf 3460sTS.Mouse.Heart.10.5.Female.fastq.gz/3460sTS.Mouse.Heart.10.5.Female_trimmed.fq.gz_stringtie.gtf 3502sTS.Mouse.Heart.12.5.Female.fastq.gz/3502sTS.Mouse.Heart.12.5.Female_trimmed.fq.gz_stringtie.gtf 5027sTS.Mouse.Heart.12.5.Female.fastq.gz/5027sTS.Mouse.Heart.12.5.Female_trimmed.fq.gz_stringtie.gtf 5028sTS.Mouse.Heart.18.5.Female.fastq.gz/5028sTS.Mouse.Heart.18.5.Female_trimmed.fq.gz_stringtie.gtf

gffcompare -R -r annotation/Mus_musculus.GRCm38.98.gtf annotation/Mus_musculus.GRCm38.98.stringtie2021.gtf -o annotation/Mus_musculus.GRCm38.98.stringtie2021.gtf.cuffcomp

python3 ensemblarize_gtf.py annotation/Mus_musculus.GRCm38.98.gtf annotation/Mus_musculus.GRCm38.98.stringtie2021.gtf.cuffcomp.annotated.gtf | sort -k1,1 -k4,4n -k5,5n > annotation/Mus_musculus.GRCm38.98.stringtie2021.fixed.gtf

mkdir annotation/Mus_musculus.GRCm38.dna.toplevel.stringtie2021.fixed 

STAR --runThreadN 4 \
--runMode genomeGenerate \
--limitGenomeGenerateRAM=141278166400 \
--genomeDir annotation/Mus_musculus.GRCm38.dna.toplevel.stringtie2021.fixed \
--genomeFastaFiles annotation/Mus_musculus.GRCm38.dna.toplevel.fa \
--sjdbGTFfile annotation/Mus_musculus.GRCm38.98.stringtie2021.fixed.gtf \
--sjdbOverhang 75

gffread annotation/Mus_musculus.GRCm38.98.stringtie2021.fixed.gtf -g annotation/Mus_musculus.GRCm38.dna.toplevel.fa -w annotation/Mus_musculus.GRCm38.98.stringtie2021.fixed.fa
gffread annotation/Mus_musculus.GRCm38.98.stringtie2021.fixed.gtf -g annotation/Mus_musculus.GRCm38.dna.toplevel.fa -x annotation/Mus_musculus.GRCm38.98.stringtie2021.fixed.CDS.fa
python3 translate_orfs.py annotation/Mus_musculus.GRCm38.98.stringtie2021.fixed.CDS.fa > annotation/Mus_musculus.GRCm38.98.stringtie2021.fixed.prot.fa

sort -k1,1 -k4,4n -k5,5n  annotation/Mus_musculus.GRCm38.98.stringtie2021.fixed.gtf > annotation/Mus_musculus.GRCm38.98.stringtie2021.fixed.sorted.gtf