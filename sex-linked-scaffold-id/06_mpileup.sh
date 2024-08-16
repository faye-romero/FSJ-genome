#!/bin/bash
#SBATCH -p rosalind
#SBATCH --job-name="mpileup"
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=20G
#SBATCH --time=5-00:00:00
#SBATCH --nodelist=bhd0056
#SBATCH --output=mpileup.log

## PART 6 SLURM script for finding sex-linked scaffolds given a female reference genome and one set of high-coverage paired-end Illumina data.##
## Faye Romero. University of Rochester. #
## 16 Dec 2023 ##

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
### Step 4. Samtools mpileup ###
################################

# Load modules
module load samtools

# Set variables
# For samtools mpileup, I want to ignore softmasked/repetitive regions by using a --positions file. Generate the --positions file using extract_unmasked_seq_pos.py

POSITIONS='/scratch/nchen11_lab/Faye/new_FSJgenome/sex_linked_scaffs/unmasked_seq_pos_FSJgenome_July2024_FINAL.bed'
REF='/scratch/nchen11_lab/Faye/new_FSJgenome/FSJgenome_July2024_FINAL_oneline.fasta'
export PATH=/scratch/nchen11_lab/Faye/new_FSJgenome:$PATH

# Create output directory, if it does not exist
mismatches_out_dir="${wd}/mismatches_filter_out"
mkdir -p "${wd}/mpileup_out"
mpileup_out_dir="${wd}/mpileup_out"

# Index bam file
cd $mismatches_out_dir
echo "$(date).....samtools index ${NAME}_aln_sorted_filt_dedup_two_mismatch.bam"

samtools index ${NAME}_aln_sorted_filt_dedup_two_mismatch.bam

cd $wd

# Execute samtools mpileup

echo "$(date).....processing ${NAME} inside ${PWD}....."

echo "$(date).....samtools mpileup --positions $POSITIONS --fasta-ref $REF ${mismatches_out_dir}/${NAME}_aln_sorted_filt_dedup_two_mismatch.bam > ${mpileup_out_dir}/${NAME}_aln_sorted_filt_dedup_two_mismatch_cov_avg.txt....."

samtools mpileup --positions $POSITIONS --fasta-ref $REF ${mismatches_out_dir}/${NAME}_aln_sorted_filt_dedup_two_mismatch.bam > ${mpileup_out_dir}/${NAME}_aln_sorted_filt_dedup_two_mismatch_cov.txt

#chgrp nchen11_lab ${mpileup_out_dir}/${NAME}_aln_sorted_filt_dedup_two_mismatch_cov.txt

echo "$(date).....done....."