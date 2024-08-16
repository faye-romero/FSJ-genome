#!/bin/bash
#SBATCH -p standard
#SBATCH --job-name="trimgalore"
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=40G
#SBATCH --time=5-00:00:00
#SBATCH --output=trimgalore.out.log
#SBATCH --error=trimgalore.err.log

## PART 1 SLURM script for finding sex-linked scaffolds given a female reference genome and one set of high-coverage paired-end Illumina data.##
## Faye Romero. University of Rochester. #
## 13 Feb 2024 ##

# Set your sample name
NAME=$(basename "$1" | cut -d'_' -f1)

# Set absolute working directory
wd=/scratch/nchen11_lab/new_FSJgenome/sex_linked_scaffs/${NAME}_out
cd $wd

# Let's begin!
echo "Processing ${NAME} within the working dir ${wd}. Let's gooooo!"
echo "$(date)."

################################
# Step 1. Trim Galore & fastqc #
################################

# Load any necessary modules
module load cutadapt/b2
module load fastqc
module load pigz
module load parallel
module load samtools
export PATH=/software/cutadapt/b2:$PATH
export PATH=/software/fastqc/0.11.8:$PATH
export PATH=/software/pigz/2.6:$PATH
export PATH=/scratch/fromero3/TrimGalore-0.6.10:$PATH
export PATH=/software/parallel/20180222:$PATH

# Set variables
TRIM='/scratch/nchen11_lab/TrimGalore-0.6.10' #path to Trimgalore folder

# Create output directory, if it does not exist
mkdir -p "${wd}/trimgalore_out"
trim_in_dir="/scratch/nchen11_lab/FSJ_WGS/rawData"
trim_out_dir="${wd}/trimgalore_out"

# Execute Trim Galore for read files
echo "$(date).....running trimgalore on $1 and $2 in ${trim_out_dir}."
echo "${TRIM}/trim_galore $1 $2 --cores 8 --quality 20 --fastqc --length 20 --paired --dont_gzip --output_dir ${trim_out_dir}"
${TRIM}/trim_galore "$1" "$2" --cores 8 --quality 20 --fastqc --length 20 --paired --dont_gzip --output_dir ${trim_out_dir}

echo "$(date).....Trim Galore has finished processing $1 and $2."

# Check if both files exist in the "trimgalore_out" subdirectory and are not empty
if [ ! -s "$trim_out_dir"/*_val_1.fq ] || [ ! -s "$trim_out_dir"/*_val_2.fq ]; then
    echo "Error: One or both of the required files (_val_1.fq and _val_2.fq) are missing or empty in $trim_out_dir directory!"
    exit 1  # Terminate the job with an error status
fi
