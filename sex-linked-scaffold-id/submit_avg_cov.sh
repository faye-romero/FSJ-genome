#!/bin/bash

# Read each sample name from the list 'list_all_samples.txt'
while IFS= read -r sample; do
    # Submit a SLURM job for each sample
    sbatch --job-name="$sample"_two_coverage --output="$sample"_two_coverage.out avg_cov.sh "$sample"
done < list_all_samples.txt

