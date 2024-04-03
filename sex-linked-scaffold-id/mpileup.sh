#!/bin/bash
#SBATCH -p rosalind
#SBATCH --job-name="mpileup"
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=30G
#SBATCH --time=5-00:00:00
#SBATCH --output=mpileup.log

## PART 6 SLURM script for finding sex-linked scaffolds given a female reference genome and one set of high-coverage paired-end Illumina data.##
## Faye Romero. University of Rochester. #
## 16 Dec 2023 ##

# To run this script, run its wrapper: sh bash_wrapper_mpileup.sh

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
### Step 4. Samtools mpileup ###
################################

# Load modules
module load samtools
# module load bedtools

# Set variables
POSITIONS='/scratch/nchen11_lab/new_FSJgenome/unmasked_seqs_FSJv3_internal_Feb2024.fasta.masked.bed' # a bed file describing the positions of unmasked (uppercase) sequences in the genome. To generate this, see extract_unmasked_seq_pos.py
REF='/scratch/nchen11_lab/new_FSJgenome/allmaps_out_v2/FSJv3_internal_Feb2024.fasta'
export PATH=/scratch/nchen11_lab/new_FSJgenome:$PATH

# Create output directory, if it does not exist
mismatches_out_dir="${wd}/mismatches_filter_out"
mkdir -p "${wd}/mpileup_out"
mpileup_out_dir="${wd}/mpileup_out"

# Execute samtools mpileup
echo "$(date).....processing ${NAME} inside ${PWD}....."

echo "$(date).....samtools mpileup --positions $POSITIONS --fasta-ref $REF ${mismatches_out_dir}/${NAME}_aln_sorted_filt_dedup_zero_mismatch.bam > ${mpileup_out_dir}/${NAME}_aln_sorted_filt_dedup_zero_mismatch_cov.txt
....."
samtools mpileup --positions $POSITIONS --fasta-ref $REF ${mismatches_out_dir}/${NAME}_aln_sorted_filt_dedup_zero_mismatch.bam > ${mpileup_out_dir}/${NAME}_aln_sorted_filt_dedup_zero_mismatch_cov.txt
echo "$(date).....done....."

echo "$(date).....samtools mpileup --positions $POSITIONS --fasta-ref $REF ${mismatches_out_dir}/${NAME}_aln_sorted_filt_dedup_two_mismatch.bam > ${mpileup_out_dir}/${NAME}_aln_sorted_filt_dedup_two_mismatch_cov.txt....."
samtools mpileup --positions $POSITIONS --fasta-ref $REF ${mismatches_out_dir}/${NAME}_aln_sorted_filt_dedup_two_mismatch.bam > ${mpileup_out_dir}/${NAME}_aln_sorted_filt_dedup_two_mismatch_cov.txt
echo "$(date).....done....."