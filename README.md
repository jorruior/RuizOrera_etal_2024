# RuizOrera_etal_2024
Set of scripts used to run the main analyses done for analyzing the evolution of genes and ORFs. **Manuscript unpublished**. We did not generate new software for this study, but used a series of published tools that are described in the folder 'Installation'.

## Installation
Any script is an standalone version and the needed dependencies are:

```
samtools==1.9
stringtie==1.2.1
price==1.0.3b
trimgalore==0.6.1
STAR==2.7.3a
gffread==0.10.1
htseq-count
bedtools==2.30.0
```
In addition, the published software ORFquant (https://github.com/lcalviell/ORFquant), PRICE (https://github.com/erhard-lab/price), and RiboseQC (https://github.com/lcalviell/Ribo-seQC) were installed to predict open reading frames.

Python3 packages
BLAST

## Data
The data generated in this study is uploaded in the European Nucleotide Archive (ENA) repository, under accession PRJEB65856.

## 1. Mapping
This folder contains a series of scripts to filter and map all generated RNA-seq and Ribo-seq data, and create the corresponding indexes. Please look at pipeline.sh to follow the different steps to map the reads to each species genome.

## 2. ORFs
This folder contains a series of scripts to predict ORFs from Ribo-seq data. ORFs were predicted using ORFquant and PRICE. Please look at pipeline.sh to follow the different steps to predict ORFs and quantify their levels of transcription and translation.

## 3. Age
This folder contains a series of scripts to predict the evolutionary age of genes and ORFs from Ribo-seq data. Please look at pipeline.sh to follow the different steps to map ORF and gene regions to other species with liftover and find homologs with BLAST.
