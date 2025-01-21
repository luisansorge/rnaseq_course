#!/bin/bash

#SBATCH --job-name="SamBamCONV"
#SBATCH --cpus-per-task=1
#SBATCH --time=1:00:00
#SBATCH --mem-per-cpu=4000M
#SBATCH --partition=pibu_el8
#SBATCH --array=0-15

# Define an array containing the sample IDs (corresponding to FASTQ files).
FASTQ_FILES=("SRR7821918" "SRR7821919" "SRR7821920" "SRR7821921" "SRR7821922" "SRR7821937" "SRR7821938" "SRR7821939" "SRR7821949" "SRR7821950" "SRR7821951" "SRR7821952" "SRR7821953" "SRR7821968" "SRR7821969" "SRR7821970")

# Select the sample ID corresponding to the current task in the job array.
XX="${FASTQ_FILES[$SLURM_ARRAY_TASK_ID]}"

# Change to directory containing the SAM files for each sample
SAM_DIR="/data/users/lansorge/rnaseq_course/mousegenome/mapped_reads"
cd ${SAM_DIR}

# Run Samtools to convert the SAM file for the current sample into a BAM file.
# Converts the input SAM file (-S) to BAM format (-b), which is a compressed binary version of the alignment.
# Directs the output to a new BAM file named after the sample ID.

apptainer exec --bind /data/ /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif samtools view -hbS ${XX}_mappedReads.sam > ${XX}_mappedReads.bam
