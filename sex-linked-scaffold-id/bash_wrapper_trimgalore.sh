#!/bin/bash

## The first wrapper script for finding sex-linked scaffolds given a female reference genome and high-coverage paired-end Illumina data for a batch of individuals. Submits one job per pair of fastq reads for each individual. ##
## Faye Romero. University of Rochester. ##
## 13 Feb 2024 ##

# Assumes that you have a parent directory with subdirectories named by SAMPLEID. Within those subdirs, you should have raw fq files labeled like so: *_1.fq.gz, *_2.fq.gz.

# This pipeline also takes in a sample list (list_all_samples.txt), which is a .txt file listing the names (SAMPLEID) of the samples you want to process, with one sample on each line. They must match the names of the subdirs mentioned above.

# Set the path to the directory containing subdirectories of raw reads per individual
data_dir="/scratch/nchen11_lab/FSJ_WGS/rawData"

# Set the desired output directory location
output_base="/scratch/nchen11_lab/new_FSJgenome/sex_linked_scaffs"

# Read the list of subdirectories from the file "list_all_samples.txt"
mapfile -t subdirectories < list_all_samples.txt

# Set the path to the SLURM script
slurm_script="trimgalore.sh"

# Loop through each subdirectory
for subdir in "${subdirectories[@]}"; do
    # Form the full path to the subdirectory
    full_subdir="${data_dir}/${subdir}"

    # Get all pairs of raw read files in the subdirectory
    read_files=("${full_subdir}"/*_1.fq.gz)  # Get all files ending with _1.fq.gz
    echo "read1 file identified: ${read_files}"
    echo "."
    echo "."
    echo "."
    for read1 in "${read_files[@]}"; do
        # Extract corresponding read2 file
        read2="${read1/_1.fq.gz/_2.fq.gz}"
		echo "read2 file identified: ${read2}"
        # Check if the corresponding read2 file exists
        if [ -e "$read2" ]; then
            # Create a unique output directory for each subdirectory in the specified location with the suffix "_out"
            output_dir="${output_base}/${subdir}_out"

            # Create the output directory only if it doesn't exist
            mkdir -p "$output_dir" && chmod 777 "$output_dir"

            # Submit a SLURM job for each pair of read files
            echo "sbatch -J ${subdir}_job -o ${output_dir}/${subdir}_trim.out.txt -e ${output_dir}/${subdir}_trim.err.txt --workdir=$output_dir $slurm_script $read1 $read2 &"
            sbatch -J "${subdir}_job" -o "${output_dir}/${subdir}_trim.out.txt" -e "${output_dir}/${subdir}_trim.err.txt" --workdir="$output_dir" "$slurm_script" "$read1" "$read2" &
        else
            echo "Error: Corresponding file not found for $read1"
        fi
    done
done

# Wait for all jobs to finish
wait
