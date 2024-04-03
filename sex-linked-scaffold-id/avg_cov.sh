#!/bin/bash
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --time=3-00:00:00
#SBATCH --mem=30G
#SBATCH --output=%x_%A.out

## PART 7 SLURM script for finding sex-linked scaffolds given a female reference genome and one set of high-coverage paired-end Illumina data.##
## Faye Romero. University of Rochester. #
## 13 Feb 2024 ##

# To run this script, run its wrapper: sh submit_avg_cov.sh

# Set variables
input_dir="/scratch/nchen11_lab/new_FSJgenome/sex_linked_scaffs"
output_dir="/scratch/nchen11_lab/new_FSJgenome/sex_linked_scaffs/avg_cov_out_two"

# Sample name passed as an argument
sample="$1"

# Navigate to the sample's directory
cd "$input_dir/${sample}_out/mpileup_out" || exit

# Calculate average coverage across each scaffold by processing the output of samtools mpileup
awk '{sum[$1]+=$4; count[$1]++} END {for (scaffold in sum) print scaffold, sum[scaffold]/count[scaffold]}' \
    "${sample}_aln_sorted_filt_dedup_two_mismatch_cov.txt" \
    > "$output_dir/${sample}_aln_sorted_filt_dedup_two_mismatch_cov_avg.txt"

