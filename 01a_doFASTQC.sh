#!/bin/bash

#SBATCH --job-name="QC"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=10:00:00
#SBATCH --mem=1000M
#SBATCH --partition=pibu_el8

# The following command runs the FastQC software using Apptainer
# which is a containerisation platform. FastQC is a tool for quality control of sequencing data.

# Runs a command within the Apptainer container, binding the `/data` directory for access.
# Specifies the container file (FastQC version 0.12.1)
# Runs FastQC, outputting results to the specified directory (/data/users/lansorge/rnaseq_course/QCoutput)
# Inputs all ".fastq.gz" files from the specified directory (/data/courses/rnaseq_course/toxoplasma_de/reads/)

apptainer exec --bind /data \
    /containers/apptainer/fastqc-0.12.1.sif \
    fastqc -o /data/users/lansorge/rnaseq_course/QCoutput /data/courses/rnaseq_course/toxoplasma_de/reads/*.fastq.gz

