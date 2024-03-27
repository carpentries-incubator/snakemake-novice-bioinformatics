# qq!/bin/bash

# This script needs the following Bioconda packages:
#   cutadapt
#   velvet
#   bbmap

#  Use cutadapt to remove adapters. Note that cutadapt works on the paired sequences
mkdir cutadapt
cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC -o cutadapt/ref_1_1.fq -p cutadapt/ref_1_2.fq reads/ref_1_1.fq reads/ref_1_2.fq
cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC -o cutadapt/ref_2_1.fq -p cutadapt/ref_2_2.fq reads/ref_2_1.fq reads/ref_2_2.fq
cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC -o cutadapt/ref_3_1.fq -p cutadapt/ref_3_2.fq reads/ref_3_1.fq reads/ref_3_2.fq

# Combine all the samples. We still want to keep the read pairs separate.
mkdir concatenate
cat cutadapt/ref_1_1.fq cutadapt/ref_2_1.fq cutadapt/ref_3_1.fq > concatenate/ref_1.fq
cat cutadapt/ref_1_2.fq cutadapt/ref_2_2.fq cutadapt/ref_3_2.fq > concatenate/ref_2.fq

# Assemble the concatenated files with Velvet. This is a 2-step process: velveth then velvetg
kmer_len=21
velveth velvet_tmp_ref_${kmer_len} ${kmer_len} -shortPaired -fastq -separate concatenate/ref_1.fq concatenate/ref_2.fq
velvetg velvet_tmp_ref_${kmer_len}
mv velvet_tmp_ref_${kmer_len}/contigs.fa contigs.fa

# Find the longest contig using the BBMap stats script
stats.sh contigs.fa | grep 'Max contig length:' > max_contig.txt
head -v max_contig.txt
