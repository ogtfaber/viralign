# viralign

Align fastq and fasta reads to virus DNA, and output a ranked score of viruses detected in the dataset.


# Pre-requirements

Installation of tools:
- bwa (https://github.com/lh3/bwa)
- samtools (https://github.com/samtools/samtools)

The scripts below assume that both tools are compiled and executable in the `./bin` directory.

## Prepare virus data

### Download virus reads

You can use any virus file, but using the one from https://rvdb.dbi.udel.edu in this example.

	mkdir -p virus
	cd ./virus
	wget https://rvdb.dbi.udel.edu/download/U-RVDBv17.0.fasta.gz
	cd ..

### Index the virus reads

	./bin/bwa index ./virus/U-RVDBv17.0.fasta.gz


# Processing the files

## Add and/or download the source data

Put the source fastq.gz files in the `./data` folder

## Run the alignment

	./bwa/bwa mem ./virus/U-RVDBv17.0.fasta.gz ./data/reads-1.fastq.gz [./data/reads-2.fastq.gz] > ./data/reads-aligned.virus.sam

## Convert to BAM

	./bin/samtools view -S -b ./data/reads-aligned.virus.sam > ./data/reads-aligned.virus.bam
	rm ./data/reads-aligned.virus.sam

## Create index

	samtools sort -o ./data/reads-aligned.virus.sorted.bam ./data/reads-aligned.virus.bam
	rm ./data/reads-aligned.virus.bam
	samtools index -b ./data/reads-aligned.virus.sorted.bam ./data/reads-aligned.virus.sorted.bam.bai

### Run the stats

	samtools idxstats ./data/reads-aligned.virus.sorted.bam > ./data/reads-aligned.virus.sorted.idxstats.csv


# Sort and filter the CSV (Python script)

	py ./bin/viralign-sort.py

## Output

Your sorted output file is at `./data/reads-aligned.virus.sorted.idxstats.viralign.csv`


# Acknowledgements

This script is developed by Onno Faber and comes without warranty of any kind. Use at your own risk. 
Initially developed for https://www.researchtothepeople.org/epithelioid-sarcoma. Thank you to all participants and organizers of these wonderful events. If you'd like to donate to *Research to the People*, visit https://www.researchtothepeople.org/donate


# Future development

I'm currently working on an environment to create a hosted version of this (and potentially other pipelines). That way nobody would have to install any environment themselves and can run open source projects like these completely wihtout any technical knowledge. If you'd like to learn more please email me at [mail@onnofaber.com](mail@onnofaber.com).

Visit https://rarematter.org to see other projects in the healthcare and rare disease space I'm working on.

