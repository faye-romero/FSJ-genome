#!/bin/bash
#SBATCH -p rosalind
#SBATCH --job-name="sambamba"
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=50G
#SBATCH --time=5-00:00:00
#SBATCH --nodelist=bhd0057
#SBATCH --output=sambamba.out.log
#SBATCH --error=sambamba.err.log

## PART 4 SLURM script for finding sex-linked scaffolds given a female reference genome and one set of high-coverage paired-end Illumina data. ##
## Faye Romero. University of Rochester. #
## 29 Nov 2023 ##

# Subdirectory path passed as an argument
SUBDIR_PATH="$1"

# Extract the sample name from the subdirectory path
NAME=$(basename "$SUBDIR_PATH" | cut -d '_' -f 1)

# Set absolute working directory
wd=/scratch/nchen11_lab/Faye/new_FSJgenome/sex_linked_scaffs/${NAME}_out
cd $wd

# Let's begin!
echo "Processing ${NAME} within the working dir ${wd}. Let's gooooo!"
echo "$(date)."


################################
## Step 4. Sambamba Markdups ###
################################

#SBATCH --job-name="sambamba.dups"
#SBATCH --output=sambamba.out.log
#SBATCH --error=sambamba.err.log

# Load conda environment and modules
module unload python python3 miniconda miniconda3
module load miniconda3
source activate sambamba

# Make sure you're in the proper working dir
cd $wd

# Set variables
sam_out_dir="${wd}/samtools_filter_out"

# Create output and temp directory, if it does not exist
mkdir -p "${wd}/sambamba_out"
sambamba_out_dir="${wd}/sambamba_out"

echo "$(date).....marking and removing dups with sambamba....."
echo "sambamba markdup ${sam_out_dir}/${NAME}_aln_sorted_filt.bam ${sambamba_out_dir}/${NAME}_aln_sorted_filt_dedup.bam \
--remove-duplicates \
--nthreads=4 \
--show-progress"

sambamba markdup ${sam_out_dir}/${NAME}_aln_sorted_filt.bam ${sambamba_out_dir}/${NAME}_aln_sorted_filt_dedup.bam \
--remove-duplicates \
--nthreads=4

#chgrp nchen11_lab ${sambamba_out_dir}/${NAME}_aln_sorted_filt_dedup.bam

echo "$(date).....done removing duplicate reads using sambamba....."
echo "."
echo "."
echo "."

# Check if file exists and is not empty
if [ ! -s "$sambamba_out_dir"/"$NAME"_aln_sorted_filt_dedup.bam ]; then
    echo "Error: The deduped file is missing or empty in $sambamba_out_dir!"
    exit 1  # Terminate the job with an error status
fi
