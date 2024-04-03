#!/bin/bash
#SBATCH -p rosalind
#SBATCH --job-name="busco.transcriptome"
#SBATCH -N 1
#SBATCH -c 12
#SBATCH --mem=20G
#SBATCH --time=5-00:00:00
#SBATCH --output=busco.transcriptome.log

## BUSCO transcriptome ##
## Faye Romero. University of Rochester. 18 March 2024 ##

# Load module, set variables, set wd
module load busco/5.2.2

BRAKER='/scratch/nchen11_lab/FSJ_GeneAnnotation/run_braker3_final'

cd /scratch/nchen11_lab/FSJ_GeneAnnotation

# Transcriptome mode: assessing assembled transcripts
echo ".....start = $(date)....."
echo "busco -m transcriptome -i ${BRAKER}/braker.codingseq -o BUSCO-auto-transcriptome --auto-lineage-euk -c 12"
busco -m transcriptome -i ${BRAKER}/braker.codingseq -o BUSCO-auto-transcriptome --auto-lineage-euk -c 12
echo ".....end = $(date)....."