#!/usr/bin/env bash

#SBATCH --job-name="SortBAM"
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=25000M
#SBATCH --time=05:00:00
#SBATCH --partition=pibu_el8
#SBATCH --array=0-15

# Define an array containing the sample IDs (corresponding to FASTQ files).
FASTQ_FILES=("SRR7821918" "SRR7821919" "SRR7821920" "SRR7821921" "SRR7821922" "SRR7821937" "SRR7821938" "SRR7821939" "SRR7821949" "SRR7821950" "SRR7821951" "SRR7821952" "SRR7821953" "SRR7821968" "SRR7821969" "SRR7821970")

# Select the sample ID corresponding to the current task in the job array.
XX="${FASTQ_FILES[$SLURM_ARRAY_TASK_ID]}"

# Change to directory containing the BAM files for each sample
BAM_DIR="/data/users/lansorge/rnaseq_course/mousegenome/mapped_reads"
cd ${BAM_DIR}

# Run Samtools to sort the BAM file for the current sample.
# Sorting organizes alignments based on genomic positions, which is necessary for downstream analyses.

# Runs the Samtools sort command.
# -m allocates 25GB of memory per thread for sorting.
# @ 4 uses 4 threads for parallel processing.
# -o specifies the name of the sorted output BAM file.
# -T temp # Specifies a prefix for temporary files created during sorting.
# Input BAM file to be sorted.

apptainer exec --bind /data/ /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif samtools sort -m 25000M -@ 4 -o ${XX}_sorted.bam -T temp ${XX}_mappedReads.bam