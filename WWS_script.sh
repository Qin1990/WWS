#!/bin/bash
#Genomic assembly
#trimmed Reads
java -jar /home/qlu/software/Trimmomatic/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 16 -trimlog tri.log WWS_1.fq.gz WWS_2.fq.gz unpaired_WWS_1.fq.gz paired_WWS_1.fq.gz unpaired_WWS_2.fq.gz paired_WWS_2.fq.gz ILLUMINACLIP:/home/qlu/software/Trimmomatic/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 MINLEN:36

#quality control
ls /home/qlu/bailachong/genome/fastqc/*gz | xargs fastqc -t 2 -o /home/qlu/bailachong/genome/fastqc/output  #paired data were used 
multiqc ./

#survey
jellyfish count -m 17 -s 20G -t 10 -o mer_counts -C WWS.fq 
jellyfish stats mer_counts
jellyfish histo -o mer_counts.histo mer_counts

#assembly
SOAPdenovo-63mer all -s contig_file -K 41 -o graph_prefix 1>ass.log 2>ass.err
GapCloser -b contig -a WWS.fasta -o WWS_genome.fasta

###########################################################################################################################################################

#Hi-C
bwa index -a bwtsw WWS_genome.fasta
bwa aln WWS_genome.fasta reads_R1.fastq.gz > reads_R1.sai
bwa aln WWS_genome.fasta reads_R2.fastq.gz > reads_R2.sai
bwa sampe WWS_genome.fasta reads_R1.sai reads_R2.sai reads_R1.fastq.gz reads_R2.fastq.gz > WWS_HiC.sam 
PreprocessSAMs.pl WWS_HiC.sam WWS_genome.fasta MBOI
filterBAM_forHiC.pl WWS_HiC.REduced.paired_only.bam WWS_HiC_clean.sam
samtools view -b WWS_HiC_clean.sam > WWS_HiC.clean.bam
cp /opt/biosoft/LACHESIS/bin/INIs/test_case.ini WWS.ini
Lachesis WWS.ini
CreateScaffoldedFasta.pl WWS_genome.fasta out/WWS
#the parameters were rewrote as following:
#species = WWS
#OUTPUT_DIR = out/90_2_3_120_10
#RE_SITE_SEQ = GCCGCGGC
#CLUSTER_N = 18
#CLUSTER_MIN_RE_SITES = 90
#CLUSTER_MAX_LINK_DENSITY = 2
#CLUSTER_NONINFORMATIVE_RATIO = 3
#ORDER_MIN_N_RES_IN_TRUNK = 120