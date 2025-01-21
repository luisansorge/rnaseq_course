#!/bin/bash

#SBATCH --job-name="MAPreads"
#SBATCH --cpus-per-task=4
#SBATCH --time=2:00:00
#SBATCH --mem-per-cpu=8000M
#SBATCH --partition=pibu_el8
#SBATCH --array=0-15

# # Define an array containing the sample IDs (corresponding to FASTQ files).
FASTQ_FILES=("SRR7821918" "SRR7821919" "SRR7821920" "SRR7821921" "SRR7821922" "SRR7821937" "SRR7821938" "SRR7821939" "SRR7821949" "SRR7821950" "SRR7821951" "SRR7821952" "SRR7821953" "SRR7821968" "SRR7821969" "SRR7821970")

# Select the sample ID corresponding to the current task in the job array.
XX="${FASTQ_FILES[$SLURM_ARRAY_TASK_ID]}"

# Define directory paths for the genome index, FASTQ input files, and output directory.
INDEX_DIR="/data/users/lansorge/rnaseq_course/mousegenome"  # Path to the directory containing the genome index
FASTQ_DIR="/data/courses/rnaseq_course/toxoplasma_de/reads" # Path to the directory containing the FASTQ files
OUTPUT_DIR=${INDEX_DIR}/mapped_reads                        # Path to the output directory for mapped reads

# Change to the genome index directory.
cd ${INDEX_DIR}

# Create the output directory if it does not already exist.
mkdir -p mapped_reads

# Run HiSAT2 using the indexed reference genome to map reads to the reference genome
# -x specifies the indexed reference genome to use for alignment.
# --rna-strandness RF Indicates the RNA sequencing data is strand-specific (reverse-forward).
# -1 specifies the input file for the first read in the pair.
# -2 specifies the input file for the second read in the pair.
# -S specifies the output SAM file for aligned reads.
# -p 4 uses 4 threads for faster processing.

apptainer exec --bind /data/ /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif hisat2 -x ${INDEX_DIR}/Mus_musculus_genome_index --rna-strandness RF -1 ${FASTQ_DIR}/${XX}_1.fastq.gz -2 ${FASTQ_DIR}/${XX}_2.fastq.gz -S ${OUTPUT_DIR}/${XX}_mappedReads.sam -p 4
 