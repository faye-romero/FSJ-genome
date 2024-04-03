#!/bin/bash
#SBATCH -p rosalind
#SBATCH --job-name="samfilter"
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=50G
#SBATCH --time=5-00:00:00
#SBATCH --output=samfilter.out.log
#SBATCH --error=samfilter.err.log
#SBATCH --mail-type=end
#SBATCH --mail-user=fromero3@ur.rochester.edu

## PART 3 SLURM script for finding sex-linked scaffolds given a female reference genome and a batch of high-coverage paired-end Illumina data from multiple individuals. ##
## Faye Romero. University of Rochester. #
## 29 Nov 2023 ##

# To run this script, run its wrapper: sh bash_wrapper_samtools_filter.sh

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
#### Step 3. Samtools filter ###
################################

# Load any necessary modules
module load samtools
module load bamtools

# Set variables

# Create output directory, if it does not exist
mkdir -p "${wd}/samtools_filter_out"
sam_out_dir="${wd}/samtools_filter_out"

# Execute samtools view (filtering)
echo "$(date).....filtering mapped reads using samtools view....."
echo "samtools view -b -f 0x2 -F 260 -q 20 -@ 4 -o ${sam_out_dir}/${NAME}_aln_sorted_temp.bam ${wd}/bwa_out/${NAME}_aln_sorted.bam"

samtools view -b -f 0x2 -F 260 -q 20 -@ 4 -o ${sam_out_dir}/${NAME}_aln_sorted_temp.bam ${wd}/bwa_out/${NAME}_aln_sorted.bam
# -b = output will be BAM
# -f 0x2 = filters for read pairs that are properly aligned (i.e., both reads in the pair mapped in the expected orientation with the correct distance between them)
# -F 260 = excludes reads that have either the 4th bit (unmapped) or the 8th bit (mate unmapped) set. This ensures that both reads in a pair are mapped.
# -q 20 = alignments with a mapping quality less than 20 are excluded.
# -@ 4 = use 4 CPU threads

echo "$(date).....sorting the filtered bam file....."
echo "samtools sort -O bam -o ${sam_out_dir}/${NAME}_aln_sorted_filt.bam ${sam_out_dir}/${NAME}_aln_sorted_temp.bam"

samtools sort -O bam -o ${sam_out_dir}/${NAME}_aln_sorted_filt.bam ${sam_out_dir}/${NAME}_aln_sorted_temp.bam

echo "$(date).....indexing the sorted, filtered bam file....."
echo "samtools index ${sam_out_dir}/${NAME}_aln_sorted_filt.bam"

cd $sam_out_dir
samtools index ${sam_out_dir}/${NAME}_aln_sorted_filt.bam
cd $wd

echo "$(date).....done filtering mapped reads and indexing the subsequent bam file....."
echo "."
echo "."
echo "."

# Check if file exists in sam_out_dir and is not empty
if [ ! -s "$sam_out_dir"/*_aln_sorted_filt.bam ] || [ ! -s "$sam_out_dir"/*_aln_sorted_filt.bam.bai ]; then
    echo "Error: _aln_sorted_filt.bam file and/or its index file is empty or doesn't exist in $sam_out_dir'!"
    exit 1  # Terminate the job with an error status
fi