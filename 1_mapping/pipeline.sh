###1. Map datasets to annotations, improve annotations and build merged datasets
#Please note that the annotation and read files should be downloaded before running this pipeline. Moreover, needed software can be found at software.txt

#Map data to Ensembl annotations
array=( "hs" "pt" "rm" "gg" "mm" "rn" )
array2=( "Homo_sapiens.GRCh38" "Pan_troglodytes.Pan_tro_3.0" "Macaca_mulatta.Mmul_10" "Gorilla_gorilla.gorGor4" "Mus_musculus.GRCm38" "Rattus_norvegicus.Rnor_6.0" )
for i in "${!array[@]}"; do
	for f in samples/Sample_*/${array[i]}*_mR*_R1*; do bash seq_pipeline.sh $f annotation/${array2[i]}.dna.toplevel annotation/${array2[i]}.98.gtf paired reverse; done
done

#Add mouse and rat developmental data from Moreira-Cardoso 2019
for f in samples/*Mouse*Heart*gz; do bash seq_pipeline.sh $f STAR/Mus_musculus.GRCm38.dna.toplevel annotation/Mus_musculus.GRCm38.98.gtf single reverse; done
for f in samples/*Rat*Heart*gz; do bash seq_pipeline.sh $f annotation/Rattus_norvegicus.Rnor_6.0.dna.toplevel annotation/Rattus_norvegicus.Rnor_6.0.98.gtf single reverse; done

#Wait until jobs are finished, then merge transcriptome
array=( "hsa" "ptr" "ggo" "mml" "mmu" "rat" )
for i in "${!array[@]}"; do qsub merge_transcriptome_${array[i]}.sh; done

#Wait until jobs are finished, then map RNA-seq and Ribo-seq to merged transcriptome
array=( "hs" "pt" "rm" "gg" "mm" "rn")
array2=( "Homo_sapiens.GRCh38" "Pan_troglodytes.Pan_tro_3.0" "Macaca_mulatta.Mmul_10" "Gorilla_gorilla.gorGor4" "Mus_musculus.GRCm38" "Rattus_norvegicus.Rnor_6.0" )
for i in "${!array[@]}"; do
	for f in samples/Sample_*/${array[i]}*_mR*_R1*; do bash seq_pipeline_sec.sh $f annotation/${array2[i]}.dna.toplevel.stringtie2021.fixed/ annotation/${array2[i]}.98.stringtie2021.fixed.gtf paired reverse; done
	for f in samples/Sample_*/${array[i]}*_mR*_R1*; do bash seq_pipeline_sec.sh $f annotation/${array2[i]}.dna.toplevel.stringtie2021.fixed/ annotation/${array2[i]}.98.stringtie2021.fixed.gtf rna29 yes; done
	for f in samples/Sample_*Ri*/${array[i]}*; do bash seq_pipeline_sec.sh $f annotation/${array2[i]}.dna.toplevel.stringtie2021.fixed/ annotation/${array2[i]}.98.stringtie2021.fixed.gtf ribo yes; done
done

#Add developmental data from Moreira-Cardoso 2019
for f in samples/*Human*gz; do bash seq_pipeline_sec.sh $f annotation/Homo_sapiens.GRCh38.dna.toplevel.stringtie2021.fixed annotation/Homo_sapiens.GRCh38.98.stringtie2021.fixed.gtf single reverse; done
for f in samples/*Macaque*gz; do bash seq_pipeline_sec.sh $f annotation/Macaca_mulatta.Mmul_10.dna.toplevel.stringtie2021.fixed annotation/Macaca_mulatta.Mmul_10.98.stringtie2021.fixed.gtf single reverse; done
for f in samples/*Mouse*gz; do bash seq_pipeline_sec.sh $f annotation/Mus_musculus.GRCm38.dna.toplevel.stringtie2021.fixed annotation/Mus_musculus.GRCm38.98.stringtie2021.fixed.gtf single reverse; done
for f in samples/Rat*gz; do bash seq_pipeline_sec.sh $f annotation/Rattus_norvegicus.Rnor_6.0.dna.toplevel.stringtie2021.fixed annotation/Rattus_norvegicus.Rnor_6.0.98.stringtie2021.fixed.gtf single reverse; done

#Get unique reads and quantify
array=( "hs" "pt" "rm" "gg" "mm" "rn")
array2=( "Homo_sapiens.GRCh38" "Pan_troglodytes.Pan_tro_3.0" "Macaca_mulatta.Mmul_10" "Gorilla_gorilla.gorGor4" "Mus_musculus.GRCm38" "Rattus_norvegicus.Rnor_6.0" )
for i in "${!array[@]}"; do
	for f in ${array[i]}*_sec/*gzAligned.sortedByCoord.out.bam; do bash extract_unique_reads.sh $f annotation/${array2[i]}.98.stringtie2021.fixed.gtf; done
	for f in ${array[i]}*_sec/*29bp*bam; do bash extract_unique_reads.sh $f annotation/${array2[i]}.98.stringtie2021.fixed.gtf; done	
done

#Wait until jobs are finished, then merge files by tissue, then merge reads (discard SB for macaque)
array=( "human" "chimp" "macaque")
array2=( "hs" "pt" "rm")
array3=( "Homo_sapiens.GRCh38" "Pan_troglodytes.Pan_tro_3.0" "Macaca_mulatta.Mmul_10" )
for i in "${!array[@]}"; do
	mkdir ${array[i]}\_lv_pooled; samtools merge -f ${array[i]}\_lv_pooled/ribo_pooled.bam ${array2[i]}\_lv_*_Ri_*_sec/*gzAligned.sortedByCoord.out.bam; samtools index ${array[i]}\_lv_pooled/ribo_pooled.bam
	mkdir ${array[i]}\_cm_pooled; samtools merge -f ${array[i]}\_cm_pooled/ribo_pooled.bam ${array2[i]}\_CM_*_Ri_*_sec/*gzAligned.sortedByCoord.out.bam; samtools index ${array[i]}\_cm_pooled/ribo_pooled.bam

	samtools merge -f ${array[i]}\_lv_pooled/rna_pooled.noxs.bam ${array2[i]}\_lv_*_mR_*_sec/*gzAligned.sortedByCoord.out.bam; samtools index ${array[i]}\_lv_pooled/rna_pooled.noxs.bam
	samtools merge -f ${array[i]}\_cm_pooled/rna_pooled.noxs.bam ${array2[i]}\_CM_*_mR_*_sec/*gzAligned.sortedByCoord.out.bam; samtools index ${array[i]}\_cm_pooled/rna_pooled.noxs.bam
	mkdir ${array[i]}\_ips_pooled; samtools merge -f ${array[i]}\_ips_pooled/rna_pooled.noxs.bam ${array2[i]}\_iPS_*_mR_*_sec/*gzAligned.sortedByCoord.out.bam; samtools index ${array[i]}\_ips_pooled/rna_pooled.noxs.bam

	samtools view -h ${array[i]}\_lv_pooled/rna_pooled.noxs.bam | awk -v strType=2 -f tagXSstrandedData.awk | samtools view -bS - > ${array[i]}\_lv_pooled/rna_pooled.bam; samtools index ${array[i]}\_cm_pooled/rna_pooled.bam
	samtools view -h ${array[i]}\_cm_pooled/rna_pooled.noxs.bam | awk -v strType=2 -f tagXSstrandedData.awk | samtools view -bS - > ${array[i]}\_cm_pooled/rna_pooled.bam; samtools index ${array[i]}\_lv_pooled/rna_pooled.bam
	samtools view -h ${array[i]}\_ips_pooled/rna_pooled.noxs.bam | awk -v strType=2 -f tagXSstrandedData.awk | samtools view -bS - > ${array[i]}\_ips_pooled/rna_pooled.bam; samtools index ${array[i]}\_ips_pooled/rna_pooled.bam

	mkdir ${array[i]}\_pooled; samtools merge -f ${array[i]}\_pooled/ribo_pooled.bam ${array[i]}\_lv_pooled/ribo_pooled.bam ${array[i]}\_cm_pooled/ribo_pooled.bam; samtools index ${array[i]}\_pooled/ribo_pooled.bam
	samtools merge -f ${array[i]}\_pooled/rna_pooled.bam ${array[i]}\_lv_pooled/rna_pooled.bam ${array[i]}\_cm_pooled/rna_pooled.bam; samtools index ${array[i]}\_pooled/rna_pooled.bam

	qsub extract_unique_reads.sh ${array[i]}\_lv_pooled/rna_pooled.bam annotation/${array3[i]}.98.stringtie2021.fixed.gtf
	qsub extract_unique_reads.sh ${array[i]}\_cm_pooled/rna_pooled.bam annotation/${array3[i]}.98.stringtie2021.fixed.gtf
	qsub extract_unique_reads.sh ${array[i]}\_ips_pooled/rna_pooled.bam annotation/${array3[i]}.98.stringtie2021.fixed.gtf

	stringtie ${array[i]}\_lv_pooled/rna_pooled.bam -e -M 0.5 -j 3 -p 4 -f 0.2 -G annotation/${array3[i]}.stringtie2021.fixed.gtf -A ${array[i]}\_lv_pooled/rna_pooled.abundance -o ${array[i]}\_lv_pooled/rna_pooled.stringtie.gtf
	stringtie ${array[i]}\_cm_pooled/rna_pooled.bam -e -M 0.5 -j 3 -p 4 -f 0.2 -G annotation/${array3[i]}.98.stringtie2021.fixed.gtf -A ${array[i]}\_cm_pooled/rna_pooled.abundance -o ${array[i]}\_cm_pooled/rna_pooled.stringtie.gtf
	stringtie ${array[i]}\_ips_pooled/rna_pooled.bam -e -M 0.5 -j 3 -p 4 -f 0.2 -G annotation/${array3[i]}.98.stringtie2021.fixed.gtf -A ${array[i]}\_ips_pooled/rna_pooled.abundance -o ${array[i]}\_ips_pooled/rna_pooled.stringtie.gtf
	stringtie ${array[i]}\_pooled/rna_pooled.bam -e -M 0.5 -j 3 -p 4 -f 0.2 -G annotation/${array3[i]}.98.stringtie2021.fixed.gtf -A ${array[i]}\_pooled/rna_pooled.abundance -o ${array[i]}\_pooled/rna_pooled.stringtie.gtf
done

##Gorilla (remove CG06)
mkdir gorilla_cm_pooled; samtools merge -f gorilla_cm_pooled/ribo_pooled.bam gg_CM_AssumboCG02___Ri_S10_R1_001.fastq.gz_sec/gg_CM_AssumboCG02___Ri_S10_R1_001_trimmed.fq.gzAligned.sortedByCoord.out.bam gg_CM_AssumboDM05___Ri_S10_R1_001.fastq.gz_sec/gg_CM_AssumboDM05___Ri_S10_R1_001_trimmed.fq.gzAligned.sortedByCoord.out.bam gg_CM_AssumboSK09___Ri_S10_R1_001.fastq.gz_sec/gg_CM_AssumboSK09___Ri_S10_R1_001_trimmed.fq.gzAligned.sortedByCoord.out.bam
samtools index gorilla_cm_pooled/ribo_pooled.bam

samtools merge -f gorilla_cm_pooled/rna_pooled.noxs.bam gg_CM_AssumboCG02___mR_S43_R1_001.fastq.gz_sec/gg_CM_AssumboCG02___mR_S43_R1_001_val_1.fq.gzAligned.sortedByCoord.out.bam gg_CM_AssumboDM05___mR_S43_R1_001.fastq.gz_sec/gg_CM_AssumboDM05___mR_S43_R1_001_val_1.fq.gzAligned.sortedByCoord.out.bam gg_CM_AssumboSK09___mR_S43_R1_001.fastq.gz_sec/gg_CM_AssumboSK09___mR_S43_R1_001_val_1.fq.gzAligned.sortedByCoord.out.bam
samtools index gorilla_cm_pooled/rna_pooled.noxs.bam
samtools merge -f gorilla_ips_pooled/rna_pooled.noxs.bam gg_iPS_Gorilla01__mR_S45_R1_001.fastq.gz_sec/gg_iPS_Gorilla01__mR_S45_R1_001_val_1.fq.gzAligned.sortedByCoord.out.bam gg_iPS_Gorilla02__mR_S45_R1_001.fastq.gz_sec/gg_iPS_Gorilla02__mR_S45_R1_001_val_1.fq.gzAligned.sortedByCoord.out.bam gg_iPS_Gorilla03__mR_S45_R1_001.fastq.gz_sec/gg_iPS_Gorilla03__mR_S45_R1_001_val_1.fq.gzAligned.sortedByCoord.out.bam
samtools index gorilla_ips_pooled/rna_pooled.noxs.bam

samtools view -h gorilla_cm_pooled/rna_pooled.noxs.bam | awk -v strType=2 -f tagXSstrandedData.awk | samtools view -bS - > gorilla_cm_pooled/rna_pooled.bam
samtools view -h gorilla_ips_pooled/rna_pooled.noxs.bam | awk -v strType=2 -f tagXSstrandedData.awk | samtools view -bS - > gorilla_ips_pooled/rna_pooled.bam
samtools index gorilla_cm_pooled/rna_pooled.bam
samtools index gorilla_ips_pooled/rna_pooled.bam

qsub extract_unique_reads.sh gorilla_cm_pooled/rna_pooled.bam annotation/Gorilla_gorilla.gorGor4.98.stringtie2021.fixed.gtf
qsub extract_unique_reads.sh gorilla_ips_pooled/rna_pooled.bam annotation/Gorilla_gorilla.gorGor4.98.stringtie2021.fixed.gtf
stringtie gorilla_cm_pooled/rna_pooled.bam -e -M 0.5 -j 3 -p 4 -f 0.2 -G annotation/Gorilla_gorilla.gorGor4.98.stringtie2021.fixed.gtf -A gorilla_cm_pooled/rna_pooled.abundance -o gorilla_cm_pooled/rna_pooled.stringtie.gtf
stringtie gorilla_ips_pooled/rna_pooled.bam -e -M 0.5 -j 3 -p 4 -f 0.2 -G annotation/Gorilla_gorilla.gorGor4.98.stringtie2021.fixed.gtf -A gorilla_ips_pooled/rna_pooled.abundance -o gorilla_ips_pooled/rna_pooled.stringtie.gtf

##Mouse
mkdir mouse_lv_pooled; samtools merge -f mouse_lv_pooled/ribo_pooled.bam mm_heart_001_Ri.fastq.gz_sec/mm_heart_001_Ri_trimmed.fq.gzAligned.sortedByCoord.out.bam mm_heart_002_Ri.fastq.gz_sec/mm_heart_002_Ri_trimmed.fq.gzAligned.sortedByCoord.out.bam mm_heart_003_Ri.fastq.gz_sec/mm_heart_003_Ri_trimmed.fq.gzAligned.sortedByCoord.out.bam mm_heart_004_Ri.fastq.gz_sec/mm_heart_004_Ri_trimmed.fq.gzAligned.sortedByCoord.out.bam mm_heart_005_Ri.fastq.gz_sec/mm_heart_005_Ri_trimmed.fq.gzAligned.sortedByCoord.out.bam mm_heart_006_Ri.fastq.gz_sec/mm_heart_006_Ri_trimmed.fq.gzAligned.sortedByCoord.out.bam
samtools index mouse_lv_pooled/ribo_pooled.bam

samtools merge -f mouse_lv_pooled/rna_pooled.noxs.bam mm_heart_001_mR_R1_raw.fastq.gz_sec/mm_heart_001_mR_R1_raw_val_1.fq.gzAligned.sortedByCoord.out.bam mm_heart_002_mR_R1_raw.fastq.gz_sec/mm_heart_002_mR_R1_raw_val_1.fq.gzAligned.sortedByCoord.out.bam mm_heart_003_mR_R1_raw.fastq.gz_sec/mm_heart_003_mR_R1_raw_val_1.fq.gzAligned.sortedByCoord.out.bam mm_heart_004_mR_R1_raw.fastq.gz_sec/mm_heart_004_mR_R1_raw_val_1.fq.gzAligned.sortedByCoord.out.bam mm_heart_005_mR_R1_raw.fastq.gz_sec/mm_heart_005_mR_R1_raw_val_1.fq.gzAligned.sortedByCoord.out.bam mm_heart_006_mR_R1_raw.fastq.gz_sec/mm_heart_006_mR_R1_raw_val_1.fq.gzAligned.sortedByCoord.out.bam
samtools index mouse_lv_pooled/rna_pooled.noxs.bam

samtools view -h mouse_lv_pooled/rna_pooled.noxs.bam | awk -v strType=2 -f tagXSstrandedData.awk | samtools view -bS - > mouse_lv_pooled/rna_pooled.bam
samtools index mouse_lv_pooled/rna_pooled.bam

qsub extract_unique_reads.sh mouse_lv_pooled/rna_pooled.bam annotation/Mus_musculus.GRCm38.98.stringtie2021.fixed.gtf
stringtie mouse_lv_pooled/rna_pooled.bam -e -M 0.5 -j 3 -p 4 -f 0.2 -G annotation/Mus_musculus.GRCm38.98.stringtie2021.fixed.gtf -A mouse_lv_pooled/rna_pooled.abundance -o mouse_lv_pooled/rna_pooled.stringtie.gtf

#Rattus_norvegicus
mkdir rat_lv_pooled; samtools merge -f rat_lv_pooled/ribo_pooled.bam rn_heart_006_Ri.fastq.gz_sec/rn_heart_006_Ri_trimmed.fq.gzAligned.sortedByCoord.out.bam rn_heart_007_Ri.fastq.gz_sec/rn_heart_007_Ri_trimmed.fq.gzAligned.sortedByCoord.out.bam rn_heart_008_Ri.fastq.gz_sec/rn_heart_008_Ri_trimmed.fq.gzAligned.sortedByCoord.out.bam rn_heart_009_Ri.fastq.gz_sec/rn_heart_009_Ri_trimmed.fq.gzAligned.sortedByCoord.out.bam rn_heart_010_Ri.fastq.gz_sec/rn_heart_010_Ri_trimmed.fq.gzAligned.sortedByCoord.out.bam
samtools index rat_lv_pooled/ribo_pooled.bam

samtools merge -f rat_lv_pooled/rna_pooled.noxs.bam rn_lv_006_ctrl_mR_S16_L002_R1_001.fastq.gz_sec/rn_lv_006_ctrl_mR_S16_L002_R1_001_val_1.fq.gzAligned.sortedByCoord.out.bam rn_lv_007_ctrl_mR_S17_L002_R1_001.fastq.gz_sec/rn_lv_007_ctrl_mR_S17_L002_R1_001_val_1.fq.gzAligned.sortedByCoord.out.bam rn_lv_008_ctrl_mR_S18_L002_R1_001.fastq.gz_sec/rn_lv_008_ctrl_mR_S18_L002_R1_001_val_1.fq.gzAligned.sortedByCoord.out.bam rn_lv_009_ctrl_mR_S19_L002_R1_001.fastq.gz_sec/rn_lv_009_ctrl_mR_S19_L002_R1_001_val_1.fq.gzAligned.sortedByCoord.out.bam rn_lv_010_ctrl_mR_S20_L002_R1_001.fastq.gz_sec/rn_lv_010_ctrl_mR_S20_L002_R1_001_val_1.fq.gzAligned.sortedByCoord.out.bam
samtools index rat_lv_pooled/rna_pooled.noxs.bam

samtools view -h rat_lv_pooled/rna_pooled.noxs.bam | awk -v strType=2 -f tagXSstrandedData.awk | samtools view -bS - > rat_lv_pooled/rna_pooled.bam
samtools index rat_lv_pooled/rna_pooled.bam

qsub extract_unique_reads.sh rat_lv_pooled/rna_pooled.bam annotation/Rattus_norvegicus.Rnor_6.0.98.stringtie2021.fixed.gtf
stringtie rat_lv_pooled/rna_pooled.bam -e -M 0.5 -j 3 -p 4 -f 0.2 -G annotation/Rattus_norvegicus.Rnor_6.0.98.stringtie2021.fixed.gtf -A rat_lv_pooled/rna_pooled.abundance -o rat_lv_pooled/rna_pooled.stringtie.gtf

#Wait until jobs are finished, filter transcriptome only with expressed transcripts
for f in ../1_mapping/*pooled/*_lvSJ.out.tab; do awk '{if ($5 < 3) print $1"\t"$2"\t"$3"\t"$7"\t"$8"\t"$4}' $f | sed 's/1$/+/' | sed 's/2$/-/' > $f.bed; done
python3 filter_exp_transcriptome.py annotation/Homo_sapiens.GRCh38.98.stringtie2021.fixed.gtf human_lv_pooled/rna_pooled.stringtie.gtf,human_cm_pooled/rna_pooled.stringtie.gtf annotation/custom/Homo_sapiens.GRCh38.98.stringtie2021.heart.gtf
python3 filter_exp_transcriptome.py annotation/Pan_troglodytes.Pan_tro_3.0.98.stringtie2021.fixed.gtf chimp_lv_pooled/rna_pooled.stringtie.gtf,chimp_cm_pooled/rna_pooled.stringtie.gtf annotation/custom/Pan_troglodytes.Pan_tro_3.0.98.stringtie2021.heart.gtf
python3 filter_exp_transcriptome.py annotation/Gorilla_gorilla.gorGor4.98.stringtie2021.fixed.gtf gorilla_cm_pooled/rna_pooled.stringtie.gtf,gorilla_cm_pooled/rna_pooled.stringtie.gtf annotation/custom/Gorilla_gorilla.gorGor4.98.stringtie2021.heart.gtf
python3 filter_exp_transcriptome.py annotation/Macaca_mulatta.Mmul_10.98.stringtie2021.fixed.gtf macaque_lv_pooled/rna_pooled.stringtie.gtf,macaque_cm_pooled/rna_pooled.stringtie.gtf annotation/custom/Macaca_mulatta.Mmul_10.98.stringtie2021.heart.gtf
python3 filter_exp_transcriptome.py annotation/Mus_musculus.GRCm38.98.stringtie2021.fixed.gtf mouse_lv_pooled/rna_pooled.stringtie.gtf,mouse_lv_pooled/rna_pooled.stringtie.gtf annotation/custom/Mus_musculus.GRCm38.98.stringtie2021.heart.gtf
python3 filter_exp_transcriptome.py annotation/Rattus_norvegicus.Rnor_6.0.98.stringtie2021.fixed.gtf rat_lv_pooled/rna_pooled.stringtie.gtf,rat_lv_pooled/rna_pooled.stringtie.gtf annotation/custom/Rattus_norvegicus.Rnor_6.0.98.stringtie2021.heart.gtf

#Build price indexes
array=( "Homo_sapiens.GRCh38" "Pan_troglodytes.Pan_tro_3.0" "Gorilla_gorilla.gorGor4" "Macaca_mulatta.Mmul_10" "Mus_musculus.GRCm38" "Rattus_norvegicus.Rnor_6.0" )
array2=( "homo_sapiens.98" "pan_troglodytes.98" "gorilla_gorilla.98" "macaca_mulatta.98" "mus_musculus.98" "rattus_norvegicus.98" )
for i in "${!array[@]}"; do
	sort -k12,12 -k3,3 annotation/custom/${array[i]}.stringtie2021.heart.gtf > annotation/custom/${array[i]}.stringtie2021.heart.sorted.gtf
	gedi -e IndexGenome -s fasta/${array[i]}.dna.toplevel.fa -a annotation/custom/${array[i]}.stringtie2021.heart.sorted.gtf -n ${array2[i]}.heart.stringtie2021 -nokallisto -nobowtie -nostar
done

