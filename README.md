# Evolution of translational control and novel open reading frames in human and non-human primate hearts.
Set of scripts used to run the main analyses done for **"Ruiz-Orera J, Miller D, Greiner J, et al. Evolution of translational control and novel open reading frames in human and non-human primate hearts" (Manuscript unpublished)**. We did not generate new software for this study, but used a series of published tools that are described in the folder 'Installation'.

## Installation
Any script is an standalone version and the needed dependencies are:

```
samtools==1.9
stringtie==2.1.1
price==1.0.3b
trimgalore==0.6.1
STAR==2.7.3a
gffread==0.10.1
htseq-count==2.0.2
bedtools==2.30.0
fastqc==0.11.9
gffcompare==0.11.6
blast==2.14.0
```
In addition, the published software ORFquant (https://github.com/lcalviell/ORFquant), PRICE (https://github.com/erhard-lab/price), and RiboseQC (https://github.com/lcalviell/Ribo-seQC) were installed to predict open reading frames. Finally, these python packages are required: Biopython (1.78), pyliftover (0.4)

## Data
The data generated in this study is uploaded in the European Nucleotide Archive (ENA) repository, under accession PRJEB65856. Other files, such as transcript and genome annotations and sequences, can be found in Ensembl (version 98 was used in this study).

## 1. Mapping
This folder contains a series of scripts to filter and map all generated RNA-seq and Ribo-seq data, and create the corresponding indexes. Please look at pipeline.sh to follow the different steps to map the reads to each species genome:
- Map RNA-seq data to Ensembl annotations.
- Use mapped information to assemble new genes and transcripts.
- Map RNA-seq and Ribo-seq data to assembled transcriptomes.
- Extract uniquely mapped reads.
- Merge mapped datasets by species and tissue.
- Filter genes and transcripts with minimum expression.
- Build PRICE indexes for ORF prediction (step 2).

## 2. ORFs
This folder contains a series of scripts to predict ORFs from Ribo-seq data. ORFs were predicted using ORFquant and PRICE. Please look at pipeline.sh to follow the different steps to predict ORFs and quantify their levels of transcription and translation:
- Extract Ribo-seq quality parameters with RiboseQC. Build ORFquant indexes.
- Run ORFquant and PRICE to predict ORFs.
- Merge ORFs per species and tissue/cell line.
- Quantify RNA-seq and Ribo-seq in predicted ORFs.

## 3. Age
This folder contains a series of scripts to predict the evolutionary age of genes and ORFs from Ribo-seq data. Please look at pipeline.sh to follow the different steps to map ORF and gene regions to other species with liftover and find homologs with BLAST:
- Build BLAST indexes.
- Find gene and transcript homologs with pyliftover.
- Find ORF homologs with pyliftover.
- Check if the counterpart ORF regions have in-frame reads in other species.
