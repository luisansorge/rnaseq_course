#!/usr/bin/env bash

#SBATCH --job-name="IndexBAM"
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4000M
#SBATCH --time=01:00:00
#SBATCH --partition=pibu_el8
#SBATCH --array=0-15

# Define an array containing the sample IDs (corresponding to FASTQ files).
FASTQ_FILES=("SRR7821918" "SRR7821919" "SRR7821920" "SRR7821921" "SRR7821922" "SRR7821937" "SRR7821938" "SRR7821939" "SRR7821949" "SRR7821950" "SRR7821951" "SRR7821952" "SRR7821953" "SRR7821968" "SRR7821969" "SRR7821970")

# Select the sample ID corresponding to the current task in the job array.
XX="${FASTQ_FILES[$SLURM_ARRAY_TASK_ID]}"

# Change to directory containing the sorted BAM files for each sample
BAM_DIR="/data/users/lansorge/rnaseq_course/mousegenome/mapped_reads"
cd ${BAM_DIR}

# Run Samtools to index the sorted BAM file for the current sample.
# Indexing creates an auxiliary .bai file for the corresponding .bam file, allowing efficient access to specific regions of the alignment data.

# Generates an index (.bai file) for the input sorted BAM file.

apptainer exec --bind /data/ /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif samtools index ${XX}_sorted.bam