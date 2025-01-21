#!/bin/bash

#SBATCH --job-name="INDEXgenome"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=5:00:00
#SBATCH --mem=8000M
#SBATCH --partition=pibu_el8

# The following command runs the Hisat2 genome indexing tool using Apptainer

# Runs a command within the Apptainer container, binding the `/data/` directory for access.
# Specifies the container file that includes Hisat2 and Samtools.
# Runs the Hisat2 genome indexing command, taking the reference genome (Mus_musculus.GRCm39.dna.primary_assembly.fa) FASTA file as input.
# Mus_musculus_genome_index specifies the base name for the generated genome index files.

apptainer exec --bind /data/ \
    /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif \
    hisat2-build /data/users/lansorge/rnaseq_course/Mus_musculus.GRCm39.dna.primary_assembly.fa Mus_musculus_genome_index
