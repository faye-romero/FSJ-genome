#!/bin/bash

## Submit script for 07_avg_cov.sh. ##
## Faye Romero. University of Rochester. ##

# Read samples from list
while IFS= read -r sample; do
    # Submit job for each sample
    sbatch --job-name="$sample"_two_coverage --output="$sample"_two_coverage.out 07_avg_cov.sh "$sample"
done < SAMPLE_LIST.txt