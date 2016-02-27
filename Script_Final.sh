#! /usr/bin/bash

#spliting the paired ends sra file
#for j in SRR311514*.sra ; do
#/scratch/evopuser/fastq-dump --split-files --gzip $j ; done

#quality of the files
#fastqc ./SRR311514*.fastq.gz

#to get genome wget "link_to _the _genome.tar.gz"
# unpack the genome file
#gunzip galGal4.fa.gz

#indexing of the genome
#STAR --runMode genomeGenerate --genomeDir star/ --genomeFastaFiles ~/hairyballs/data/genome/galGal4.fa

#mapping
for i in *_1.fastq.gz ; do
STAR --genomeDir ./star --readFilesIn $i ${i/_1/_2} --outFileNamePrefix ${i%_1*}.star. --outSAMattributes All --runThreadN 8 \
--readFilesCommand zcat; done


#from sam to bam and sort file

#for k in *.star.Aligned.out.sam ; do
#samtools view -bS $k | samtools sort -o $k.bam; done


#remove sam files
#rm *.star.Aligned.out.sam

#index bam file
#for j in *.star.Aligned.out.sam.bam; do
#samtools index $j; done

#create table with counts
/scratch/evopuser/tools/subread-1.4.0-p1/bin/featureCounts -t exon -g gene_id -a ~/hairyballs/data/genome/galGal4.gtf -o RefSeq.counts140 SRR3115140.star.Aligned.out.bam






