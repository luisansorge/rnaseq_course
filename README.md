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

_Table 1: Sample names and corresponding experimental conditions. Sequencing was performed in paired-end mode, resulting in two FASTQ files per sample, representing read 1 and read 2, respectively._

## 1. Quality Control
**Quality control assessment of raw reads using FastQC:**

_run script 01a_doFASTQC.sh_

**Creation of FastQC reports for all samples, enabling the possibility for MultiQC analysis to compare the quality of reads across all samples:**

_run script 01b_doMULTIQC.sh_

## 2. Mapping reads to the reference genome
**Indexing of the [Mus musculus reference genome](https://ftp.ensembl.org/pub/release-113/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz) using HISAT2:**

_run script 02_doINDEX.sh_ 

**Mapping raw reads to the indexed reference genome using HISAT2:**

_run script 03_doMAPreads.sh_

**Conversion of HISAT2 output SAM files to BAM files using Samtools:**

_run script 04_doSAMBAMconversion.sh_

**Sorting BAM files to organise alignments by genomic position using Samtools, necessary for downstream analyses:**

_run script 05_doBAMsort.sh_

**Indexing sorted BAM files with Samtools to create auxiliary .bai files for efficient access to specific regions of alignment data:**

_run script 06_doBAMindex.sh_

**Following this workflow produces an indexed alignment of each sample to the Mus musculus reference genome.**

## 3. Counting the number of reads per gene
**Counting the number of reads assigned to exons to enable comparison of read counts across different genes, providing a foundation for subsequent differential expression analysis using featureCounts:**

_run script 07_doCOUNTreads.sh_

**Creation of a summary report of featureCounts results for comparing the assignment of reads across all samples with MultiQC:**

_run script 07b_doMULTIQCfeaturecounts.sh_

## 4. Exploratory, Differential Expression, and Overrepresentation Analysis
**Locally performed on RStudio, principal component analysis (PCA) and differential gene expression analysis with pairwise comparisons were conducted using the DESeq2 package. Additionally, Gene Ontology (GO) terms containing more differentially expressed genes than expected were identified using clusterProfiler:**

_run script 08_DESeq2_DE_analysis in RStudio_
