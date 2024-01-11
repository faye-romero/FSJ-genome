#!/bin/bash
#SBATCH -p rosalind
#SBATCH --job-name="cutadapt"
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH --time=13-00:00:00
#SBATCH --output=cutadapt.log
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=fromero3@ur.rochester.edu

## Script to run Cutadapt on multiple sets of raw PacBio Hifi reads ##
## Faye Romero, University of Rochester ##
## 2 May 2022 ##

# Set wd
cd /scratch/nchen11_lab/FSJgenome/HiFi

# Load modules
module load cutadapt/b2

# Run Cutadapt on raw read sets
for SAMPLE in 083029 190552 224304
do
echo "date Run cut adapt for ${SAMPLE}"
cutadapt --discard-trimmed --overlap=35 -b ATCTCTCTCAACAACAACAACGGAGGAGGAGGAAAAGAGAGAGAT -b ATCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGAT -o ./hifi_cutadapt.${SAMPLE}.fastq.gz /scratch/nchen11_lab/pacbio/hifi_reads.${SAMPLE}.fastq.gz
done

# Combine raw reads into one fastq file
# cat hifi_cutadapt.*.fastq.gz > FSJ.hifi.cut.fastq.gz