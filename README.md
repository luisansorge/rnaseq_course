# Detection of differentially expressed genes from bulk RNA-sequencing data
This workflow begins with Illumina sequencing data in the form of FASTQ files and aims to generate lists of differentially expressed (DE) genes between two experimental groups. Additionally, it seeks to identify gene ontology (GO) terms that are significantly enriched among the DE genes.

The dataset comprises 3-5 biological replicates from two distinct mice tissues, comparing Toxoplasma-infected specimens with uninfected controls. These samples represent a subset of the data published by [Singhania et al. (2019)](https://pmc.ncbi.nlm.nih.gov/articles/PMC6599044/), with raw sequencing files (FASTQ format) obtained from the Gene Expression Omnibus (GEO) repository (accession number: GSE119855). The RNA-seq libraries were prepared using a strand-specific protocol and subsequently sequenced on an Illumina HiSeq 4000 platform in paired-end configuration.

## Sample list

|Sample	|Group|
|-------|------|
|SRR7821921|	Lung Infected|
SRR7821922|	Lung Infected
SRR7821918	|Lung Infected
SRR7821919	|Lung Infected
SRR7821920|	Lung Infected
SRR7821937	|Lung Control
SRR7821938	|Lung Control
SRR7821939	|Lung Control
SRR7821949	|Blood Infected
SRR7821950	|Blood Infected
SRR7821951	|Blood Infected
SRR7821952	|Blood Infected
SRR7821953	|Blood Infected
SRR7821968	|Blood Control
SRR7821969	|Blood Control
SRR7821970|	Blood Control

_Table 1: shows the sample names and corresponding experimental conditions. The sequencing was performed in paired-end mode, resulting in two FASTQ files per sample, representing read 1 and read 2, respectively._

## 1. Quality Control
**Quality control assessment using the raw reads using fastqc:**

_run script 01a_doFASTQC.sh_

**Creation of fastqc reports for all samples, enabling the possibility for MultiQC analysis which creates a report of the fastqc results in order to compare the quality of reads across all samples:**

_run script 01b_doMULTIQC.sh_

## 2. Mapping reads to the reference genome
**Indexing of the [Mus musculus reference genome](https://ftp.ensembl.org/pub/release-113/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz) using HISAT2:**

_run script 02_doINDEX.sh_ 

**Mapping of the raw reads to the indexed reference genome using HISAT2:**

_run script 03_doMAPreads.sh_

**HISAT2 output in the format of SAM files, conversion to BAM files through Samtools:**

_run script 04_doSAMBAMconversion.sh_

**Sorting of the BAM files in order to organise alignments based on genomic position with Samtools, which is necessary for downstream analyses:**

_run script 05_doBAMsort.sh_

**Indexing the sorted BAM files with Samtools, creating an auxiliary .bai file for the corresponding .bam file, allowing efficient access to specific regions of the alignment data:**

_run script 06_doBAMindex.sh_

**Following this workflow allows for the production of an indexed alignement of each sample to the reference genome of Mus musculus.**

## 3. Counting the number of reads per gene
**Counting the numbers of reads that assigned to exons enabling the comparison of read counts across different genes and providing a foundation for subsequent differential expression analysis using featureCounts:**

_run script 07_doCOUNTreads.sh_

**Creation of a summary report of the featureCounts results in order to compare the assignment of reads across all samples with MultiQC:**

_run script 07b_doMULTIQCfeaturecounts.sh_
