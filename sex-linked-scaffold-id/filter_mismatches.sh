#!/bin/bash
#SBATCH -p rosalind
#SBATCH --job-name="filter.mismatches"
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=25G
#SBATCH --time=5-00:00:00
#SBATCH --output=filter.mismatches.log

## PART 5 SLURM script for finding sex-linked scaffolds given a female reference genome and a batch of high-coverage paired-end Illumina data from multiple individuals. ##
## Faye Romero. University of Rochester. #
## 9 Dec 2023 ##

# To run this script, run its wrapper: sh bash_wrapper_filter_mismatches.sh

# Subdirectory path passed as an argument
SUBDIR_PATH="$1"

# Extract the sample name from the subdirectory path
NAME=$(basename "$SUBDIR_PATH" | cut -d '_' -f 1)

# Set absolute working directory
wd=/scratch/nchen11_lab/new_FSJgenome/sex_linked_scaffs/${NAME}_out
cd $wd

# Let's begin!
echo "Processing ${NAME} within the working dir ${wd}. Let's gooooo!"
echo "$(date)."

################################
#### Step 3. Bamtools filter ###
################################

# Load any necessary modules
module load samtools
module load bamtools

# Set variables

# Create output directory, if it does not exist
sambamba_out_dir="${wd}/sambamba_out"
mkdir -p "${wd}/mismatches_filter_out"
mis_out_dir="${wd}/mismatches_filter_out"

# Execute bamtools

# Zero mismatches allowed (strict filtering)
echo "$(date).....filtering for zero mismatches (strict) using bamtools....."
echo "bamtools filter -in ${sambamba_out_dir}/${NAME}_aln_sorted_filt_dedup.bam -tag NM:<=0 -out ${mis_out_dir}/${NAME}_aln_sorted_filt_dedup_zero_mismatch.bam"

bamtools filter -in ${sambamba_out_dir}/${NAME}_aln_sorted_filt_dedup.bam -tag "NM:<=0" -out ${mis_out_dir}/${NAME}_aln_sorted_filt_dedup_zero_mismatch.bam

echo "$(date).....finished filtering for zero mismatches (strict) using bamtools....."
echo "."
echo "."
echo "."

# Two mismatches allowed (intermediate filtering)
echo "$(date).....filtering for two mismatches (intermediate) using bamtools....."
echo "bamtools filter -in ${sambamba_out_dir}/${NAME}_aln_sorted_filt_dedup.bam -tag NM:<=2 -out ${mis_out_dir}/${NAME}_aln_sorted_filt_dedup_two_mismatch.bam"

bamtools filter -in ${sambamba_out_dir}/${NAME}_aln_sorted_filt_dedup.bam -tag "NM:<=2" -out ${mis_out_dir}/${NAME}_aln_sorted_filt_dedup_two_mismatch.bam

echo "$(date).....finished filtering for two mismatches (intermediate) using bamtools....."
echo "."
echo "."
echo "."
