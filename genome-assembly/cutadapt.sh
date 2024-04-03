#!/bin/bash
#SBATCH --partition=rosalind --time=13-00:00:00 --output=cutadapt.2May2022.log --mail-type=all
#SBATCH --mem=64G -c 8

## Cutadapt ##
## Nancy Chen, Faye Romero. University of Rochester. 02 May 2022 ##

# Navigate to working directory
cd /scratch/nchen11_lab/FSJgenome/HiFi

# Load cutadapt
module load cutadapt/b2

# Cutadapt removes reads with adapters

for SAMPLE in 083029 190552 224304
do
echo "date Run cut adapt for ${SAMPLE}"
cutadapt --discard-trimmed --overlap=35 -b ATCTCTCTCAACAACAACAACGGAGGAGGAGGAAAAGAGAGAGAT -b ATCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGAT -o ./hifi_cutadapt.${SAMPLE}.fastq.gz /scratch/nchen11_lab/pacbio/hifi_reads.${SAMPLE}.fastq.gz
done