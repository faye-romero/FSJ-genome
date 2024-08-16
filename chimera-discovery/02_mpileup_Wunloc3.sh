#!/bin/bash
#SBATCH -p rosalind
#SBATCH --job-name="mpileup_Wunloc3"
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=20G
#SBATCH --nodelist=bhd0057
#SBATCH --time=1-00:00:00
#SBATCH --output=mpileup_Wunloc3.log

## Calculate coverage for all positions across a single scaffold (in my case, the chimeric scaffold ChrW_unloc_scaf_3) from a filtered bam file ##
## Faye Romero. University of Rochester. #
## 15 July 2024 ##

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
REF='/scratch/nchen11_lab/Faye/new_FSJgenome/FSJgenome_July2024_FINAL.fasta'
export PATH=/scratch/nchen11_lab/Faye/new_FSJgenome:$PATH

# Create output directory, if it does not exist
mismatches_out_dir="${wd}/mismatches_filter_out"
mkdir -p "${wd}/mpileup_out"
mpileup_out_dir="${wd}/mpileup_out"

# Index bam file
# cd $mismatches_out_dir
# echo "$(date).....samtools index ${NAME}_aln_sorted_filt_dedup_two_mismatch.bam"

# samtools index ${NAME}_aln_sorted_filt_dedup_two_mismatch.bam

# cd $wd

# Execute samtools mpileup for all positions across ChrW_unloc_scaf_3
echo "$(date).....processing ${NAME} inside ${PWD}....."

echo "$(date).....samtools mpileup -a -r ChrW_unloc_scaf_3 --fasta-ref $REF ${mismatches_out_dir}/${NAME}_aln_sorted_filt_dedup_two_mismatch.bam > /scratch/nchen11_lab/Faye/new_FSJgenome/sex_linked_scaffs/Wunloc3_cov_AllPos/${NAME}_aln_sorted_filt_dedup_two_mismatch_cov_Wunloc3_AllPos.txt....."

samtools mpileup -a -r ChrW_unloc_scaf_3 --fasta-ref $REF ${mismatches_out_dir}/${NAME}_aln_sorted_filt_dedup_two_mismatch.bam > /scratch/nchen11_lab/Faye/new_FSJgenome/sex_linked_scaffs/Wunloc3_cov_AllPos/${NAME}_aln_sorted_filt_dedup_two_mismatch_cov_Wunloc3_AllPos.txt
echo "$(date).....done....."