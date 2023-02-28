#!/bin/bash
#$ -cwd -R yes -l m_mem_free=60G -l h_rt=60:00:00 -pe smp 4
#$ -o logs/mapping2.o -e logs/mapping2.e

source ~/.bash_profile


fastq=$1 #Files should end in '.fastq.gz'
star=$2
gtf=$3
class=$4 #single,paired,ribo
orient=$5 

if [ $class == "single" ]
then
  bf1="$(basename -- $fastq)"
  ff1="${bf1//\.fastq.gz/_trimmed.fq.gz}"
  echo "Filtering and trimming single RNA-seq "$bf1
  trim_galore --length 25 -q 30 --trim-n --gzip $fastq -o $bf1\_sec

  fastqc $bf1\_sec/$ff1
  echo "Mapping "$ff1
  STAR --runThreadN 4 \
  --genomeDir $star \
  --twopassMode Basic \
  --alignSJDBoverhangMin 6 \
  --alignSJoverhangMin 500 \
  --readFilesIn \
  $bf1\_sec/$ff1 \
  --readFilesCommand zcat \
  --outFilterMismatchNmax 2 \
  --outFilterMultimapNmax 5 \
  --outSAMattributes All \
  --outSAMtype BAM SortedByCoordinate \
  --quantMode GeneCounts \
  --outFileNamePrefix $bf1\_sec/$ff1 \
  --limitOutSJcollapsed 10000000 \
  --limitIObufferSize=300000000 \
  --outSAMattrIHstart 0 \
  --outSAMstrandField intronMotif \
  --outFilterIntronMotifs RemoveNoncanonical \  
  --outFilterType BySJout \
  --outTmpKeep none

  samtools index $bf1\_sec/$ff1\Aligned.sortedByCoord.out.bam

  echo "Counting, reverse by default"
  htseq-count -f bam -s $orient $bf1\_sec/$ff1\Aligned.sortedByCoord.out.bam $gtf > $bf1\_sec/$ff1.htseq.unique.counts
  htseq-count -f bam -s $orient -t CDS $bf1\_sec/$ff1\Aligned.sortedByCoord.out.bam $gtf > $bf1\_sec/$ff1.htseq.cds.unique.counts

elif [ $class == "paired" ]
then
  fastq1=$1
  fastq2=${fastq1/_R1/_R2}

  bf1="$(basename -- $fastq1)"
  bf2="$(basename -- $fastq2)"
  ff1="${bf1//.fastq.gz/_val_1.fq.gz}"
  ff2="${bf2//.fastq.gz/_val_2.fq.gz}"

  echo "Filtering and trimming paired RNA-seq "$fastq1" "$fastq2
  trim_galore --paired --length 25 -q 30 --gzip $fastq1 $fastq2 -o $bf1\_sec

  fastqc $bf1\_sec/$ff1
  fastqc $bf1\_sec/$ff2
  echo "Mapping "$ff1" "$ff2
  STAR --runThreadN 4 \
  --genomeDir $star \
  --twopassMode Basic \
  --alignSJDBoverhangMin 6 \
  --alignSJoverhangMin 500 \
  --readFilesIn \
  $bf1\_sec/$ff1 \
  $bf1\_sec/$ff2 \
  --readFilesCommand zcat \
  --outFilterMismatchNmax 4 \
  --outFilterMultimapNmax 20 \
  --outSAMattributes All \
  --outSAMtype BAM SortedByCoordinate \
  --quantMode GeneCounts \
  --outFileNamePrefix $bf1\_sec/$ff1 \
  --limitOutSJcollapsed 10000000 \
  --limitIObufferSize=300000000 \
  --outSAMattrIHstart 0 \
  --outSAMstrandField intronMotif \
  --outFilterIntronMotifs RemoveNoncanonical \
  --outFilterType BySJout \
  --outTmpKeep none

  samtools index $bf1\_sec/$ff1\Aligned.sortedByCoord.out.bam

  echo "Counting, reverse by default"
  htseq-count -f bam -r pos -s $orient $bf1\_sec/$ff1\Aligned.sortedByCoord.out.bam $gtf > $bf1\_sec/$ff1.htseq.unique.counts
  htseq-count -f bam -r pos -s $orient -t CDS $bf1\_sec/$ff1\Aligned.sortedByCoord.out.bam $gtf > $bf1\_sec/$ff1.htseq.cds.unique.counts

elif [ $class == "ribo" ]
then
  bf1="$(basename -- $fastq)"
  ff1="${bf1//.fastq.gz/_trimmed.fq.gz}"
  echo "Filtering and trimming single Ribo-seq "$bf1
  trim_galore --length 25 --trim-n --gzip $fastq -o $bf1\_sec

  echo "Filtering contaminants "$ff1
  bowtie2 --seedlen=25 -p 4 --un $bf1\_sec/$ff1\_filtered.fastq -x contaminants $bf1\_sec/$ff1>/dev/null

  fastqc $bf1\_sec/$ff1
  echo "Mapping "$ff1
  STAR --runThreadN 4 \
  --genomeDir $star \
  --twopassMode Basic \
  --alignSJDBoverhangMin 6 \
  --alignSJoverhangMin 500 \
  --readFilesIn \
  $bf1\_sec/$ff1\_filtered.fastq \
  --outFilterMismatchNmax 2 \
  --outFilterMultimapNmax 5 \
  --outSAMattributes All \
  --outSAMtype BAM SortedByCoordinate \
  --quantMode GeneCounts \
  --outFileNamePrefix $bf1\_sec/$ff1 \
  --limitOutSJcollapsed 10000000 \
  --limitIObufferSize=300000000 \
  --outFilterType BySJout \
  --outTmpKeep none

  STAR --runThreadN 4 \
  --genomeDir $star \
  --twopassMode Basic \
  --readFilesIn \
  $bf1\_sec/$ff1\_filtered.fastq \
  --outFilterMismatchNmax 2 \
  --outFilterMultimapNmax 1 \
  --outSAMattributes All \
  --outSAMtype BAM SortedByCoordinate \
  --quantMode GeneCounts \
  --outFileNamePrefix $bf1\_sec/$ff1\_forprice \
  --limitOutSJcollapsed 10000000 \
  --limitIObufferSize=300000000 \
  --outFilterType BySJout \
  --alignSJoverhangMin 1000 \
  --outTmpKeep none \
  --alignEndsType EndToEnd

  samtools index $bf1\_sec/$ff1\_forpriceAligned.sortedByCoord.out.bam
  samtools index $bf1\_sec/$ff1\Aligned.sortedByCoord.out.bam
  echo "Stringtie"
  stringtie $bf1\_sec/$ff1\Aligned.sortedByCoord.out.bam -M 0.5 -j 3 -p 4 -f 0.2 -G $3 -e -A $bf1\_sec/$ff1\_abundance -o $bf1\_sec/$ff1\_stringtie.gtf

  echo "Counting, reverse by default"
  htseq-count -f bam -s $orient $bf1\_sec/$ff1\Aligned.sortedByCoord.out.bam $gtf > $bf1\_sec/$ff1.htseq.unique.counts
  htseq-count -f bam -s $orient -t CDS $bf1\_sec/$ff1\Aligned.sortedByCoord.out.bam $gtf > $bf1\_sec/$ff1.htseq.cds.unique.counts

elif [ $class == "rna29" ]
then
  bf1="$(basename -- $fastq)"
  ff1="${bf1//.fastq.gz/_trimmed.fq.gz}"
  ff2="${ff1//\.fq.gz/.29bp_5prime.fq.gz}"
  echo "Filtering and trimming single RNA-seq 29bp "$bf1
  trim_galore --length 29 --trim-n --gzip $fastq -o $bf1\_29bp_sec
  trim_galore --hardtrim5 29 --gzip $bf1\_29bp_sec/$ff1 -o $bf1\_29bp_sec

  echo "Filtering contaminants "$ff2
  bowtie2 --seedlen=25 -p 4 --un $bf1\_29bp_sec/$ff2\_filtered.fastq -x contaminants $bf1\_29bp_sec/$ff2>/dev/null

  echo "Mapping "$ff2
  STAR --runThreadN 4 \
  --genomeDir $star \
  --twopassMode Basic \
  --alignSJDBoverhangMin 6 \
  --alignSJoverhangMin 500 \
  --readFilesIn \
  $bf1\_29bp_sec/$ff2\_filtered.fastq \
  --outFilterMismatchNmax 2 \
  --outFilterMultimapNmax 5 \
  --outSAMattributes All \
  --outSAMtype BAM SortedByCoordinate \
  --quantMode GeneCounts \
  --outFileNamePrefix $bf1\_29bp_sec/$ff2 \
  --limitOutSJcollapsed 10000000 \
  --limitIObufferSize=300000000 \
  --outFilterType BySJout \
  --outTmpKeep none  

  samtools index $bf1\_29bp_sec/$ff2\Aligned.sortedByCoord.out.bam
  echo "Stringtie"
  stringtie $bf1\_29bp_sec/$ff2\Aligned.sortedByCoord.out.bam -M 0.5 -j 3 -p 4 -f 0.2 -G $3 -e -A $bf1\_29bp_sec/$ff2\_abundance -o $bf1\_29bp_sec/$ff2\_stringtie.gtf

  echo "Counting, reverse by default"
  htseq-count -f bam -s $orient $bf1\_29bp_sec/$ff2\Aligned.sortedByCoord.out.bam $gtf > $bf1\_29bp_sec/$ff2.htseq.unique.counts
  htseq-count -f bam -s $orient -t CDS $bf1\_29bp_sec/$ff2\Aligned.sortedByCoord.out.bam $gtf > $bf1\_29bp_sec/$ff2.htseq.cds.unique.counts
fi

