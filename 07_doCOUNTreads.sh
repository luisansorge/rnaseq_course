#!/usr/bin/env bash

#SBATCH --job-name="COUNTreads"
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=1000M
#SBATCH --time=05:00:00
#SBATCH --partition=pibu_el8

# Change to directory containing the sorted BAM files for each sample
BAM_DIR="/data/users/lansorge/rnaseq_course/mousegenome/mapped_reads"
cd ${BAM_DIR}

# Define the output directory for storing the gene counts
OUTPUT_DIR="/data/users/lansorge/rnaseq_course/gene_counts"

# Define the directory containing the annotation file (GTF format) for the reference genome.
ANNOTATION_DIR="/data/users/lansorge/rnaseq_course/mousegenome"

# Run featureCounts to quantify gene expression based on the mapped reads.
# -p enables paired-end read counting by treating read pairs as a single unit.
# -s 2 specifies reverse strandedness.
# -Q 10 sets the minimum mapping quality of reads to 10.
# -t exon specifies the feature type to count in this case exons.
# -g gene_id specifies the attribute in the GTF file to use for grouping in this case the gene_id.
# -a supplies the annotation file (GTF) for the reference genome.
# -o specifies the output file for the combined gene counts.
# *_sorted.bam inputs all sorted BAM files in the directory for read counting.

apptainer exec --bind /data/ /containers/apptainer/subread_2.0.1--hed695b0_0.sif featureCounts -p -s 2 -Q 10 -t exon -g gene_id -a ${ANNOTATION_DIR}/Mus_musculus.GRCm39.113.gtf -o ${OUTPUT_DIR}/combined_counts_quality.txt *_sorted.bam