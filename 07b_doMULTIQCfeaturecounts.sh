#!/bin/bash

#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=1000M
#SBATCH --time=01:00:00
#SBATCH --partition=pibu_el8

# Load MultiQC to create a report of the featureCounts results in order to compare the assignment of reads across all samples
module load MultiQC/1.11-foss-2021a

# Inputs all files from the featureCounts analysis in the specified directory (/data/users/lansorge/rnaseq_course/gene_counts)
FEATURECOUNTS_RESULTS_DIR="/data/users/lansorge/rnaseq_course/gene_counts"
MULTIQC_OUTPUT_DIR="/data/users/lansorge/rnaseq_course/multiQC_output"

# The following command runs MultiQC
# -o: defines the output directory
# --interactive enables the creation of a HTML report

multiqc $FEATURECOUNTS_RESULTS_DIR -o $MULTIQC_OUTPUT_DIR --interactive --title "MultiQC Report - featureCounts"