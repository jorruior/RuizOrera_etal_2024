#Build BLAST transcript and ORF indexes
for f in "Homo_sapiens.GRCh38.98" "Pan_troglodytes.Pan_tro_3.0.98" "Gorilla_gorilla.gorGor4.98" "Macaca_mulatta.Mmul_10.98" "Mus_musculus.GRCm38.98" "Rattus_norvegicus.Rnor_6.0.98"; do 
	makeblastdb -dbtype nucl -in ../1_mapping/annotation/$f.stringtie2021.fixed.fa -out ../1_mapping/annotation/$f.stringtie2021.fixed
done
for f in ../1_mapping/*_pooled/*_plus_NTG.nucl.fa; do makeblastdb -dbtype nucl -in $f -out $f.db; done
for f in ../1_mapping/*_pooled/*.orfs.fa; do makeblastdb -dbtype prot -in $f -out $f.db; done
intersectBed -wo -s -b <(grep -P "\tCDS\t" ../1_mapping/annotation/Homo_sapiens.GRCh38.98.gtf) -a ../1_mapping/human_pooled/human.orfs.bed | awk '{print $4"\t"$26}' | sort -u | sed 's/"//g' | sed 's/;//' | awk '{a[$1]=a[$1] FS $2} END{for(i in a) print i a[i]}' | sed 's/ /#/' | sed 's/ /;/g' | sed 's/#/\t/' > ../1_mapping/annotation/human_ensembl98_alloverlaps.txt

#Find ORFs in the transcriptomes of other species
gffread <(sed 's/transcript_id/trans_id/' ../1_mapping/human_pooled/human.orfs.gtf | sed 's/orf_id/transcript_id/') -g /fast/AG_Huebner/Jorge/GENOMES/fasta/Homo_sapiens.GRCh38.dna.toplevel.fa -x ../1_mapping/human_pooled/human.orfs.cds.fa

#1. Find transcript homologs
bash find_gene_homologs.sh ../annotation/Homo_sapiens.GRCh38.98.stringtie2021.fixed.gtf ../annotation/Homo_sapiens.GRCh38.98.gtf ../annotation/Homo_sapiens.GRCh38.98.stringtie2021.fixed.fa human
bash find_gene_homologs.sh ../annotation/Pan_troglodytes.Pan_tro_3.0.98.stringtie2021.fixed.gtf ../annotation/Pan_troglodytes.Pan_tro_3.0.98.gtf ../annotation/Pan_troglodytes.Pan_tro_3.0.98.stringtie2021.fixed.fa chimp
bash find_gene_homologs.sh ../annotation/Gorilla_gorilla.gorGor4.98.stringtie2021.fixed.gtf ../annotation/Gorilla_gorilla.gorGor4.98.gtf ../annotation/Gorilla_gorilla.gorGor4.98.stringtie2021.fixed.fa gorilla
bash find_gene_homologs.sh ../annotation/Macaca_mulatta.Mmul_10.98.stringtie2021.fixed.gtf ../annotation/Macaca_mulatta.Mmul_10.98.gtf ../annotation/Macaca_mulatta.Mmul_10.98.stringtie2021.fixed.fa macaque
bash find_gene_homologs.sh ../annotation/Mus_musculus.GRCm38.98.stringtie2021.fixed.gtf ../annotation/Mus_musculus.GRCm38.98.gtf ../annotation/Mus_musculus.GRCm38.98.stringtie2021.fixed.fa mouse
bash find_gene_homologs.sh ../annotation/Rattus_norvegicus.Rnor_6.0.98.stringtie2021.fixed.gtf ../annotation/Rattus_norvegicus.Rnor_6.0.98.gtf ../annotation/Rattus_norvegicus.Rnor_6.0.98.stringtie2021.fixed.fa rat

cat liftover/human_transcriptome_to_*.exons.status.txt > liftover/human_transcriptome.exons.status.txt
cat liftover/chimp_transcriptome_to_*.exons.status.txt > liftover/chimp_transcriptome.exons.status.txt
cat liftover/macaque_transcriptome_to_*.exons.status.txt > liftover/macaque_transcriptome.exons.status.txt
cat liftover/mouse_transcriptome_to_*.exons.status.txt > liftover/mouse_transcriptome.exons.status.txt
cat liftover/rat_transcriptome_to_*.exons.status.txt > liftover/rat_transcriptome.exons.status.txt


#2. Find ORF homologs
for f in "human" "chimp" "macaque"; do 
	bash find_orf_homologs.sh ../1_mapping/$f\_pooled/$f.orfs.gtf_Psites.bed ../1_mapping/$f\_pooled/$f.orfs.fa $f
done
for f in "mouse" "rat"; do 
	bash find_orf_homologs.sh ../1_mapping/$f\_lv_pooled/$f\_lv.orfs.gtf_Psites.bed ../1_mapping/$f\_lv_pooled/$f\_lv.orfs.fa $f\_lv
done
bash find_orf_homologs.sh ../1_mapping/gorilla_cm_pooled/gorilla_cm.orfs.gtf_Psites.bed ../1_mapping/gorilla_cm_pooled/gorilla_cm.orfs.fa gorilla_cm
cat liftover/human_*synteny.txt > liftover/human.orfs_to_transcripts_heart.synteny.txt
cat liftover/chimp_*synteny.txt > liftover/chimp.orfs_to_transcripts_heart.synteny.txt
cat liftover/gorilla_*synteny.txt > liftover/gorilla.orfs_to_transcripts_heart.synteny.txt
cat liftover/macaque_*synteny.txt > liftover/macaque.orfs_to_transcripts_heart.synteny.txt
cat liftover/mouse_*synteny.txt > liftover/mouse.orfs_to_transcripts_heart.synteny.txt
cat liftover/rat_*synteny.txt > liftover/rat.orfs_to_transcripts_heart.synteny.txt


#3. Evaluate and quantify in-frame aligned regions in the counterpart regions of ORFs in other species.
for f in "human" "chimp" "macaque"; do 
	bash find_region_homologs.sh ../1_mapping/$f\_pooled/$f.orfs.bed ../1_mapping/$f\_pooled/$f.orfs.gtf_Psites.bed $f
done
