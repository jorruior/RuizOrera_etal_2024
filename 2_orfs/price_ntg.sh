#!/bin/bash
#$ -cwd -R yes -l m_mem_free=200G -l h_rt=96:00:00 -pe smp 1
#$ -o logs/price.o -e logs/price.e

bam=$1
species=$2 #homo_sapiens.98 or mus_musculus.98
exp=$3
expg=$4

gedi -e Price -reads $bam -prefix $bam -genomic $species -progress
awk '{if ($9 < 0.05) print $3"\t"$2"-"$5"\t"$1}' $bam.orfs.tsv | sed 's/:/\t/' | sed 's/+\t/\t+\t/' | sed 's/-\t/\t-\t/' | awk '{split($3,a,"|"); $3="";for (i in a){print a[i]"\t"$0}}' | sed 's/-/\t/' | awk '{if (length($6) > 1) print $3"\tPrice\tCDS\t"$1+1"\t"$2"\t.\t"$4"\t.\tORF_id \""$5"\"; gene_id \""$6"\";"}' > $bam.price.gtf

if [ $species == "homo_sapiens.98.heart.stringtie2021" ]
then
    python3 add_nearcognate_price.py $bam.price.gtf ../1_mapping/annotation/Homo_sapiens.GRCh38.98.stringtie2021.fixed.fa ../1_mapping/fasta/Homo_sapiens.GRCh38.dna.toplevel.fa ../1_mapping/annotation/Homo_sapiens.GRCh38.98.stringtie2021.fixed.prot.fa ../1_mapping/annotation/Homo_sapiens.GRCh38.98.stringtie2021.fixed.gtf ../1_mapping/annotation/Homo_sapiens.GRCh38.98.gtf_psites.bed $exp $expg
fi

if [ $species == "pan_troglodytes.98.heart.stringtie2021" ]
then
    python3 add_nearcognate_price.py $bam.price.gtf ../1_mapping/annotation/Pan_troglodytes.Pan_tro_3.0.98.stringtie2021.fixed.fa ../1_mapping/fasta/Pan_troglodytes.Pan_tro_3.0.dna.toplevel.fa ../1_mapping/annotation/Pan_troglodytes.Pan_tro_3.0.98.stringtie2021.fixed.prot.fa ../1_mapping/annotation/Pan_troglodytes.Pan_tro_3.0.98.stringtie2021.fixed.gtf ../1_mapping/annotation/Pan_troglodytes.Pan_tro_3.0.98.gtf_psites.bed $exp $expg
fi

if [ $species == "gorilla_gorilla.98.heart.stringtie2021" ]
then
    python3 add_nearcognate_price.py $bam.price.gtf ../1_mapping/annotation/Gorilla_gorilla.gorGor4.98.stringtie2021.fixed.fa ../1_mapping/fasta/Gorilla_gorilla.gorGor4.dna.toplevel.fa ../1_mapping/annotation/Gorilla_gorilla.gorGor4.98.stringtie2021.fixed.prot.fa ../1_mapping/annotation/Gorilla_gorilla.gorGor4.98.stringtie2021.fixed.gtf ../1_mapping/annotation/Gorilla_gorilla.gorGor4.98.gtf_psites.bed $exp $expg
fi

if [ $species == "macaca_mulatta.98.heart.stringtie2021" ]
then
    python3 add_nearcognate_price.py $bam.price.gtf ../1_mapping/annotation/Macaca_mulatta.Mmul_10.98.stringtie2021.fixed.fa ../1_mapping/fasta/Macaca_mulatta.Mmul_10.dna.toplevel.fa ../1_mapping/annotation/Macaca_mulatta.Mmul_10.98.stringtie2021.fixed.prot.fa ../1_mapping/annotation/Macaca_mulatta.Mmul_10.98.stringtie2021.fixed.gtf ../1_mapping/annotation/Macaca_mulatta.Mmul_10.98.gtf_psites.bed $exp $expg
fi

if [ $species == "mus_musculus.98.heart.stringtie2021" ]
then
    guixa
    python3 add_nearcognate_price.py $bam.price.gtf ../1_mapping/annotation/Mus_musculus.GRCm38.98.stringtie2021.fixed.fa ../1_mapping/fasta/Mus_musculus.GRCm38.dna.toplevel.fa ../1_mapping/annotation/Mus_musculus.GRCm38.98.stringtie2021.fixed.prot.fa ../1_mapping/annotation/Mus_musculus.GRCm38.98.stringtie2021.fixed.gtf ../1_mapping/annotation/Mus_musculus.GRCm38.98.gtf_psites.bed $exp $expg
fi

if [ $species == "rattus_norvegicus.98.heart.stringtie2021" ]
then
    guixa
    python3 add_nearcognate_price.py $bam.price.gtf ../1_mapping/annotation/Rattus_norvegicus.Rnor_6.0.dna.toplevel.stringtie2021.fixed.fa ../1_mapping/fasta/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa ../1_mapping/annotation/Rattus_norvegicus.Rnor_6.0.dna.toplevel.stringtie2021.fixed.prot.fa ../1_mapping/annotation/Rattus_norvegicus.Rnor_6.0.98.stringtie2021.fixed.gtf ../1_mapping/annotation/Rattus_norvegicus.Rnor_6.0.98.gtf_psites.bed $exp $expg
fi

if [ $species == "heterocephalus_glaber.98.heart.stringtie2021" ]
then
    guixa
    python3 add_nearcognate_price.py $bam.price.gtf ../1_mapping/annotation/Heterocephalus_glaber_male.HetGla_1.0.dna.toplevel.stringtie2021.fixed.fa ../1_mapping/fasta/Heterocephalus_glaber_male.HetGla_1.0.dna.toplevel.fa ../1_mapping/annotation/Heterocephalus_glaber_male.HetGla_1.0.dna.toplevel.stringtie2021.fixed.prot.fa ../1_mapping/annotation/Heterocephalus_glaber_male.HetGla_1.0.98.stringtie2021.fixed.gtf ../1_mapping/annotation/Heterocephalus_glaber_male.HetGla_1.0.98.gtf_psites.bed $exp $expg
fi

if [ $species == "homo_sapiens.98.blt.stringtie2020" ]
then
    python3 add_nearcognate_price.py $bam.price.gtf ../1_mapping/annotation/Homo_sapiens.GRCh38.98.stringtie2020.fixed.fa ../1_mapping/fasta/Homo_sapiens.GRCh38.dna.toplevel.fa ../1_mapping/annotation/Homo_sapiens.GRCh38.98.stringtie2020.fixed.prot.fa ../1_mapping/annotation/Homo_sapiens.GRCh38.98.stringtie2020.fixed.gtf ../1_mapping/annotation/Homo_sapiens.GRCh38.98.gtf_psites.bed $exp $expg
fi

if [ $species == "macaca_mulatta.98.blt.stringtie2020" ]
then
    python3 add_nearcognate_price.py $bam.price.gtf ../1_mapping/annotation/Macaca_mulatta.Mmul_10.98.stringtie2020.fixed.fa ../1_mapping/fasta/Macaca_mulatta.Mmul_10.dna.toplevel.fa ../1_mapping/annotation/Macaca_mulatta.Mmul_10.98.stringtie2020.fixed.prot.fa ../1_mapping/annotation/Macaca_mulatta.Mmul_10.98.stringtie2020.fixed.gtf ../1_mapping/annotation/Macaca_mulatta.Mmul_10.98.gtf_psites.bed $exp $expg
fi

if [ $species == "mus_musculus.98.blt.stringtie2020" ]
then
    guixa
    python3 add_nearcognate_price.py $bam.price.gtf ../1_mapping/annotation/Mus_musculus.GRCm38.98.stringtie2020.fixed.fa ../1_mapping/fasta/Mus_musculus.GRCm38.dna.toplevel.fa ../1_mapping/annotation/Mus_musculus.GRCm38.98.stringtie2020.fixed.gtf.CDS.prot.fa ../1_mapping/annotation/Mus_musculus.GRCm38.98.stringtie2020.fixed.gtf ../1_mapping/annotation/Mus_musculus.GRCm38.98.gtf_psites.bed $exp $expg
fi

if [ $species == "gallus_gallus.98.blt.stringtie2020" ]
then
    guixa
    python3 add_nearcognate_price.py $bam.price.gtf ../1_mapping/annotation/Gallus_gallus.GRCg6a.98.stringtie2020.fixed.fa ../1_mapping/fasta/Gallus_gallus.GRCg6a.dna.toplevel.fa ../1_mapping/annotation/Gallus_gallus.GRCg6a.98.stringtie2020.fixed.prot.fa ../1_mapping/annotation/Gallus_gallus.GRCg6a.98.stringtie2020.fixed.gtf ../1_mapping/annotation/Gorilla_gorilla.gorGor4.98.gtf_psites.bed $exp $expg
fi

if [ $species == "homo_sapiens.98" ]
then
    python3 add_nearcognate_price.py $bam.price.gtf ../1_mapping/annotation/Homo_sapiens.GRCh38.cdna_and_ncrna.ENST.all.fa ../1_mapping/fasta/Homo_sapiens.GRCh38.dna.toplevel.fa ../1_mapping/annotation/Homo_sapiens.GRCh38.pep.ENST.all.fa ../1_mapping/annotation/Homo_sapiens.GRCh38.98_forORFquant.gtf ../1_mapping/annotation/Homo_sapiens.GRCh38.98.gtf_psites.bed $exp $expg
fi
