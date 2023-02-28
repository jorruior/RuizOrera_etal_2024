#Build indexes for ORFquant
array=( "hsa" "ptr" "ggo" "mml" "mmu" "rat" )
array2=( "Homo_sapiens.GRCh38.98" "Pan_troglodytes.Pan_tro_3.0.98" "Gorilla_gorilla.gorGor4.98" "Macaca_mulatta.Mmul_10.98" "Mus_musculus.GRCm38.98" "Rattus_norvegicus.Rnor_6.0.98" )
array3=( "hsa38.98" "Pantro3.0.98" "gorGor4.98" "Mmul10.98" "mmu.98" "rat.98" )
for i in "${!array[@]}"; do
	Rscript riboseqc_primates.R X ${array[i]} 1 0 ../1_mapping/annotation/${array2[i]}.stringtie2021.fixed.gtf ${array3[i]}.stringtie2021 X
done

for i in "${!array[@]}"; do
	Rscript riboseqc_primates.R X ${array[i]} 1 0 ../1_mapping/annotation/custom/${array2[i]}.stringtie2021.heart.gtf ${array3[i]}.heart.stringtie2021 X
done

#Wait until they are finished, then run ORFquant
for f in ../1_mapping/hs*_Ri_*sec/*gzAligned.sortedByCoord.out.bam; do qsub -N hs orfquant_primates.sh $f hsa ../1_mapping/custom/Homo_sapiens.GRCh38.98.stringtie2021.heart.gtf hsa38.98.heart.stringtie2021 ../1_mapping/custom/Homo_sapiens.GRCh38.98.stringtie2021.heart.gtf_Rannot none none; done
for f in ../1_mapping/pt*_Ri_*sec/*gzAligned.sortedByCoord.out.bam; do qsub -N pt orfquant_primates.sh $f ptr ../1_mapping/custom/Pan_troglodytes.Pan_tro_3.0.98.stringtie2021.heart.gtf Pantro3.0.98.heart.stringtie2021 ../1_mapping/custom/Pan_troglodytes.Pan_tro_3.0.98.stringtie2021.heart.gtf_Rannot none none; done
for f in ../1_mapping/gg*_Ri_*sec/*gzAligned.sortedByCoord.out.bam; do qsub -N gg orfquant_primates.sh $f ggo ../1_mapping/custom/Gorilla_gorilla.gorGor4.98.stringtie2021.heart.gtf gorGor4.98.heart.stringtie2021 ../1_mapping/custom/Gorilla_gorilla.gorGor4.98.stringtie2021.heart.gtf_Rannot none none; done
for f in ../1_mapping/rm*_Ri_*sec/*gzAligned.sortedByCoord.out.bam; do qsub -N rm orfquant_primates.sh $f mml ../1_mapping/custom/Macaca_mulatta.Mmul_10.98.stringtie2021.heart.gtf Mmul10.98.heart.stringtie2021 ../1_mapping/custom/Macaca_mulatta.Mmul_10.98.stringtie2021.heart.gtf_Rannot none none; done
for f in ../1_mapping/mm*_Ri*sec/*gzAligned.sortedByCoord.out.bam; do qsub -N mm orfquant_primates.sh $f mmu ../1_mapping/custom/Mus_musculus.GRCm38.98.stringtie2021.heart.gtf mmu.98.heart.stringtie2021 ../1_mapping/custom/Mus_musculus.GRCm38.98.stringtie2021.heart.gtf_Rannot none none; done
for f in ../1_mapping/rn*_Ri*sec/*gzAligned.sortedByCoord.out.bam; do qsub -N rn orfquant_primates.sh $f rat ../1_mapping/custom/Rattus_norvegicus.Rnor_6.0.98.stringtie2021.heart.gtf rat.98.heart.stringtie2021 ../1_mapping/custom/Rattus_norvegicus.Rnor_6.0.98.stringtie2021.heart.gtf_Rannot none none; done

#Combined orfquant
qsub -pe smp 7 orfquant_primates.sh ../1_mapping/human_pooled/ribo_pooled.bam hsa ../1_mapping/custom/Homo_sapiens.GRCh38.98.stringtie2021.heart.gtf hsa38.98.heart.stringtie2021 ../1_mapping/custom/Homo_sapiens.GRCh38.98.stringtie2021.heart.gtf_Rannot ../1_mapping/human_pooled/rna_pooled.stringtie.gtf ../1_mapping/human_pooled/rna_pooled.abundance
qsub -pe smp 7 orfquant_primates.sh ../1_mapping/chimp_pooled/ribo_pooled.bam ptr ../1_mapping/custom/Pan_troglodytes.Pan_tro_3.0.98.stringtie2021.heart.gtf Pantro3.0.98.heart.stringtie2021 ../1_mapping/custom/Pan_troglodytes.Pan_tro_3.0.98.stringtie2021.heart.gtf_Rannot ../1_mapping/chimp_pooled/rna_pooled.stringtie.gtf ../1_mapping/chimp_pooled/rna_pooled.abundance
qsub -pe smp 7 orfquant_primates.sh ../1_mapping/macaque_pooled/ribo_pooled.bam mml ../1_mapping/custom/Macaca_mulatta.Mmul_10.98.stringtie2021.heart.gtf Mmul10.98.heart.stringtie2021 ../1_mapping/custom/Macaca_mulatta.Mmul_10.98.stringtie2021.heart.gtf_Rannot ../1_mapping/macaque_pooled/rna_pooled.stringtie.gtf ../1_mapping/macaque_pooled/rna_pooled.abundance

array=( "hs" "pt" "gg" "rm" "mm" "rn" )
array2=( "homo_sapiens.98.heart.stringtie2021" "pan_troglodytes.98.heart.stringtie2021" "gorilla_gorilla.98.heart.stringtie2021" "macaca_mulatta.98.heart.stringtie2021" "mus_musculus.98.heart.stringtie2021" "rattus_norvegicus.98.heart.stringtie2021" )
for i in "${!array[@]}"; do
	for f in ../1_mapping/${array[i]}*_Ri*sec/*forpriceAligned.sortedByCoord.out.bam; do qsub price_ntg.sh $f ${array2[i]} none none; done
done

#Wait until they are finished, next merge ORFs and create unified lists of ORFs
array=( "human_cm" "human_lv" "human" "chimp_cm" "chimp_lv" "chimp" "gorilla_cm" "macaque_cm" "macaque_lv" "macaque" "mouse_lv" "rat_lv" )
array2=( "Homo_sapiens.GRCh38.98" "Homo_sapiens.GRCh38.98" "Homo_sapiens.GRCh38.98" "Pan_troglodytes.Pan_tro_3.0.98" "Pan_troglodytes.Pan_tro_3.0.98" "Pan_troglodytes.Pan_tro_3.0.98" "Gorilla_gorilla.gorGor4.98" "Macaca_mulatta.Mmul_10.98" "Macaca_mulatta.Mmul_10.98" "Macaca_mulatta.Mmul_10.98" "Mus_musculus.GRCm38.98" "Rattus_norvegicus.Rnor_6.0.98" )
for i in "${!array[@]}"; do
	qsub ORF_mapper_plus_consensus.sh ${array[i]} ${array2[i]}
done
#Wait until they are finished

#Quantify RNA-seq in ORFs
array=( "hs" "pt" "gg" "rm" "mm" "rn" )
array2=( "human" "chimp" "gorilla_cm" "macaque" "mouse_lv" "rat_lv" )
cd ../1_mapping
for i in "${!array[@]}"; do
	for f in ${array[i]}_*mR_*sec/*29bp*bam; do qsub htseq_count_orf.sh $f ${array2[i]}\_pooled/${array2[i]}.orfs.gtf yes ${array2[i]}; done
done