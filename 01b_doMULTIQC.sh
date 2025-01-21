#!/bin/bash

#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=1000M
#SBATCH --time=01:00:00
#SBATCH --partition=pibu_el8

# Load MultiQC to create a report of the fastqc results in order to compare the quality of reads across all samples
module load MultiQC/1.11-foss-2021a

# Inputs all files from the fastq analysis in the specified directory (/data/users/lansorge/rnaseq_course/QCoutput)
FASTQC_RESULTS_DIR="/data/users/lansorge/rnaseq_course/QCoutput"
MULTIQC_OUTPUT_DIR="/data/users/lansorge/rnaseq_course/multiQC_output"

#Create a directory to store the multiQC report
mkdir $MULTIQC_OUTPUT_DIR

# The following command runs MultiQC
# -o: defines the output directory
# --interactive enables the creation of a HTML report

multiqc $FASTQC_RESULTS_DIR -o $MULTIQC_OUTPUT_DIR --interactive --title "MultiQC Report - Raw Data"