#!/bin/bash
#$ -cwd -R yes -l m_mem_free=120G -l h_rt=96:00:00 -pe smp 3
#$ -o logs/orfquant.o -e logs/orfquant.e

bam=$1
species=$2
gtf=$3
annotname=$4
annot=$5
exp=$6
expg=$7

Rscript riboseqc_primates.R $bam $species 7 -1 $gtf $annotname $annot
Rscript riboseqc_primates.R $bam $species 7 -2 $gtf $annotname $annot

if [ $species == "hsa" ]
then
    guixa
    python3 fix_fasta_orfquant.py $bam\_ORFquant_Detected_ORFs.gtf $bam\_ORFquant_Protein_sequences.fasta ../1_mapping/annotation/Homo_sapiens.GRCh38.98.stringtie2021.fixed.fa ../1_mapping/fasta/Homo_sapiens.GRCh38.dna.toplevel.fa ../1_mapping/annotation/Homo_sapiens.GRCh38.98.stringtie2021.fixed.prot.fa ../1_mapping/annotation/Homo_sapiens.GRCh38.98.stringtie2021.fixed.gtf ../1_mapping/annotation/Homo_sapiens.GRCh38.98.gtf_psites.bed $exp $expg
fi

if [ $species == "ptr" ]
then
    guixa
    python3 fix_fasta_orfquant.py $bam\_ORFquant_Detected_ORFs.gtf $bam\_ORFquant_Protein_sequences.fasta ../1_mapping/annotation/Pan_troglodytes.Pan_tro_3.0.98.stringtie2021.fixed.fa ../1_mapping/fasta/Pan_troglodytes.Pan_tro_3.0.dna.toplevel.fa ../1_mapping/annotation/Pan_troglodytes.Pan_tro_3.0.98.stringtie2021.fixed.prot.fa ../1_mapping/annotation/Pan_troglodytes.Pan_tro_3.0.98.stringtie2021.fixed.gtf ../1_mapping/annotation/Pan_troglodytes.Pan_tro_3.0.98.gtf_psites.bed $exp $expg
fi

if [ $species == "ggo" ]
then
    guixa
    python3 fix_fasta_orfquant.py $bam\_ORFquant_Detected_ORFs.gtf $bam\_ORFquant_Protein_sequences.fasta ../1_mapping/annotation/Gorilla_gorilla.gorGor4.98.stringtie2021.fixed.fa ../1_mapping/fasta/Gorilla_gorilla.gorGor4.dna.toplevel.fa ../1_mapping/annotation/Gorilla_gorilla.gorGor4.98.stringtie2021.fixed.prot.fa ../1_mapping/annotation/Gorilla_gorilla.gorGor4.98.stringtie2021.fixed.gtf ../1_mapping/annotation/Gorilla_gorilla.gorGor4.98.gtf_psites.bed $exp $expg
fi

if [ $species == "mml" ]
then
    guixa
    python3 fix_fasta_orfquant.py $bam\_ORFquant_Detected_ORFs.gtf $bam\_ORFquant_Protein_sequences.fasta ../1_mapping/annotation/Macaca_mulatta.Mmul_10.98.stringtie2021.fixed.fa ../1_mapping/fasta/Macaca_mulatta.Mmul_10.dna.toplevel.fa ../1_mapping/annotation/Macaca_mulatta.Mmul_10.98.stringtie2021.fixed.prot.fa ../1_mapping/annotation/Macaca_mulatta.Mmul_10.98.stringtie2021.fixed.gtf ../1_mapping/annotation/Macaca_mulatta.Mmul_10.98.gtf_psites.bed $exp $expg
fi

if [ $species == "mmu" ]
then
    guixa
    python3 fix_fasta_orfquant.py $bam\_ORFquant_Detected_ORFs.gtf $bam\_ORFquant_Protein_sequences.fasta ../1_mapping/annotation/Mus_musculus.GRCm38.98.stringtie2021.fixed.fa ../1_mapping/fasta/Mus_musculus.GRCm38.dna.toplevel.fa ../1_mapping/annotation/Mus_musculus.GRCm38.98.stringtie2021.fixed.prot.fa ../1_mapping/annotation/Mus_musculus.GRCm38.98.stringtie2021.fixed.gtf ../1_mapping/annotation/Mus_musculus.GRCm38.98.gtf_psites.bed $exp $expg
fi

if [ $species == "rat" ]
then
    guixa
    python3 fix_fasta_orfquant.py $bam\_ORFquant_Detected_ORFs.gtf $bam\_ORFquant_Protein_sequences.fasta ../1_mapping/annotation/Rattus_norvegicus.Rnor_6.0.dna.toplevel.stringtie2021.fixed.fa ../1_mapping/fasta/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa ../1_mapping/annotation/Rattus_norvegicus.Rnor_6.0.dna.toplevel.stringtie2021.fixed.prot.fa ../1_mapping/annotation/Rattus_norvegicus.Rnor_6.0.98.stringtie2021.fixed.gtf ../1_mapping/annotation/Rattus_norvegicus.Rnor_6.0.98.gtf_psites.bed $exp $expg
fi

if [ $species == "het" ]
then
    guixa
    python3 fix_fasta_orfquant.py $bam\_ORFquant_Detected_ORFs.gtf $bam\_ORFquant_Protein_sequences.fasta ../1_mapping/annotation/Heterocephalus_glaber_male.HetGla_1.0.dna.toplevel.stringtie2021.fixed.fa ../1_mapping/fasta/Heterocephalus_glaber_male.HetGla_1.0.dna.toplevel.fa ../1_mapping/annotation/Heterocephalus_glaber_male.HetGla_1.0.dna.toplevel.stringtie2021.fixed.prot.fa ../1_mapping/annotation/Heterocephalus_glaber_male.HetGla_1.0.98.stringtie2021.fixed.gtf ../1_mapping/annotation/Heterocephalus_glaber_male.HetGla_1.0.98.gtf_psites.bed $exp $expg
fi

if [ $species == "gal" ]
then
    guixa
    python3 fix_fasta_orfquant.py $bam\_ORFquant_Detected_ORFs.gtf $bam\_ORFquant_Protein_sequences.fasta ../1_mapping/annotation/Gallus_gallus.GRCg6a.98.stringtie2020.fixed.fa ../1_mapping/fasta/Gallus_gallus.GRCg6a.dna.toplevel.fa ../1_mapping/annotation/Gallus_gallus.GRCg6a.98.stringtie2020.fixed.prot.fa ../1_mapping/annotation/Gallus_gallus.GRCg6a.98.stringtie2020.fixed.gtf ../1_mapping/annotation/Gallus_gallus.GRCg6a.98.gtf_psites.bed $exp $expg
fi

if [ $species == "hsa98" ]
then
    guixa
    python3 fix_fasta_orfquant.py $bam\_ORFquant_Detected_ORFs.gtf $bam\_ORFquant_Protein_sequences.fasta ../1_mapping/annotation/Homo_sapiens.GRCh38.cdna_and_ncrna.ENST.all.fa ../1_mapping/fasta/Homo_sapiens.GRCh38.dna.toplevel.fa ../1_mapping/annotation/Homo_sapiens.GRCh38.pep.ENST.all.fa ../1_mapping/annotation/Homo_sapiens.GRCh38.98_forORFquant.gtf ../1_mapping/annotation/Homo_sapiens.GRCh38.98.gtf_psites.bed $exp $expg
fi

