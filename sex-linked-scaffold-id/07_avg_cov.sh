#!/bin/bash
#SBATCH --partition=rosalind
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --time=3-00:00:00
#SBATCH --mem=20G
#SBATCH --nodelist=bhd0056
#SBATCH --output=%x_%A.out

## PART 7 SLURM script for finding sex-linked scaffolds given a female reference genome and one set of high-coverage paired-end Illumina data.##
## Faye Romero. University of Rochester. #
## 16 Dec 2023 ##

# Set paths to input parent dir and output dir
input_dir="/scratch/nchen11_lab/Faye/please/sex_linked_scaffs"
output_dir="/scratch/nchen11_lab/Faye/please/sex_linked_scaffs/avg_cov"

# Sample name passed as an argument
sample="$1"

# Navigate to the sample's mpileup directory
cd "$input_dir/${sample}_out/mpileup_out" || exit

# Iterate through mpileup coverage file and calculate the average cov per scaffold
awk '{sum[$1]+=$4; count[$1]++} END {for (scaffold in sum) print scaffold, sum[scaffold]/count[scaffold]}' \
    "${sample}_aln_sorted_filt_dedup_two_mismatch_cov.txt" \
    > "$output_dir/${sample}_aln_sorted_filt_dedup_two_mismatch_cov_avg.txt"

