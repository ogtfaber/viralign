#!/bin/bash

# Run the alignment

~/bin/bwa/bwa mem /ref/U-RVDBv17.0.fasta.gz /data/sample.fasta > /data/reads-aligned.virus.sam

# Convert to BAM

~/bin/samtools/samtools view -S -b /data/reads-aligned.virus.sam > /data/reads-aligned.virus.bam
rm /data/reads-aligned.virus.sam

# Create index

~/bin/samtools/samtools sort -o /data/reads-aligned.virus.sorted.bam /data/reads-aligned.virus.bam
rm /data/reads-aligned.virus.bam
~/bin/samtools/samtools index -b /data/reads-aligned.virus.sorted.bam /data/reads-aligned.virus.sorted.bam.bai

# Run the stats

~/bin/samtools/samtools idxstats /data/reads-aligned.virus.sorted.bam > /data/reads-aligned.virus.sorted.idxstats.csv

head -1000 /data/reads-aligned.virus.sorted.idxstats.csv

# Test
# head -1000 /data/sample.fasta
