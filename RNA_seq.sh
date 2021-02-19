#!/bin/bash
#trimmed Reads
java -jar /home/qlu/software/Trimmomatic/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 16 -trimlog tri.log FF1_1.fq.gz FF1_2.fq.gz unpaired_FF1_1.fq.gz paired_FF1_1.fq.gz unpaired_FF1_2.fq.gz paired_FF1_2.fq.gz ILLUMINACLIP:/home/qlu/software/Trimmomatic/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 MINLEN:36
java -jar /home/qlu/software/Trimmomatic/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 16 -trimlog tri.log FF2_1.fq.gz FF2_2.fq.gz unpaired_FF2_1.fq.gz paired_FF2_1.fq.gz unpaired_FF2_2.fq.gz paired_FF2_2.fq.gz ILLUMINACLIP:/home/qlu/software/Trimmomatic/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 MINLEN:36
java -jar /home/qlu/software/Trimmomatic/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 16 -trimlog tri.log FF3_1.fq.gz FF3_2.fq.gz unpaired_FF3_1.fq.gz paired_FF3_1.fq.gz unpaired_FF3_2.fq.gz paired_FF3_2.fq.gz ILLUMINACLIP:/home/qlu/software/Trimmomatic/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 MINLEN:36
#all raw data (SF, FM, SM) processed with the same parameters

#quality control
ls /home/qlu/bailachong/genome/fastqc/*gz | xargs fastqc -t 2 -o /home/qlu/bailachong/genome/fastqc/output  #paired data were used 
multiqc ./

#align to genome
hisat2-build –p 4 WWS_genome.fa genome
hisat2 -p 16 -x ./qlu/genome -1 FF1_1.fastq,FF2_1.fastq,FF3_1.fastq,SF1_1.fastq,SF2_1.fastq,SF3_1.fastq,FM1_1.fastq,FM2_1.fastq,FM3_1.fastq,SM1_1.fastq,SM2_1.fastq,SM3_1.fastq -2 FF1_2.fastq,FF2_2.fastq,FF3_2.fastq,SF1_2.fastq,SF2_2.fastq,SF3_2.fastq,FM1_2.fastq,FM2_2.fastq,FM3_2.fastq,SM1_2.fastq,SM2_2.fastq,SM3_2.fastq –S WWS.sam
samtools view -b WWS.sam > WWS.bam
samtools sort WWS.bam -o WWS_sort.bam
samtools index WWS_sorted.bam

#identify counts
GTF="/home/qlu/bailachong/RNA_seq/hisat2/bailachong.gtf"
featureCounts -T 5 -t exon -g gene_id -a $GTF -p -o WWS_RNA.txt *sort.bam

#identify DEGs
options(stringsAsFactors=F)
#the data was imported to excel for identify the DEGs
library("DESeq2")
database <- read.table(file = "WWS_RNA.xls", sep = "\t", header = T, row.names = 1)
countData <- database[,1:12]
condition <- factor(c(rep("FF",3),rep("SF",3),rep("FM",3),rep("SM",3)), levels = c("FF", "SF","FM","SM"))
countData <- round(as.matrix(countData))
coldata <- data.frame(row.names = colnames(countData), condition)
dds <- DESeqDataSetFromMatrix(countData, DataFrame(condition), design= ~ condition )
dds <- dds[ rowSums(counts(dds)) > 1, ]
dds <- DESeq(dds)
res0vs1 <- results(dds, contrast = c("condition","FF","SF"))
table(res0vs1$pvalue <0.01)
res0vs1 <- res0vs1[order(res0vs1$pvalue),]
write.table(as.data.frame(res0vs1), file="res0vs1.xls", sep="\t", quote = F)
#other groups (FF VS FM, FF VS SM, SF VS FM, SF VS SM, FM VS SM) were processed with same parameters
