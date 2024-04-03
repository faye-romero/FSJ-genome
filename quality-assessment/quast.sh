#!/bin/bash
#SBATCH -p rosalind
#SBATCH --job-name="quast"
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=50G
#SBATCH --time=1-00:00:00
#SBATCH --output=quast.09252023.log

## QUAST ##
## Faye Romero. University of Rochester. 25 September 2023 ##

# Load modules
module load python3/3.7.1
module load quast/5.0.2
module load gcc/9.1.0
module load gnuplot/5.2.0b2
module load gnutls/3.5.15
module load zlib/1.2.11

# Set variables
DRAFT='AphCoe_V3_internal_Mar2024.fasta'
TOPDIR='/scratch/nchen11_lab/new_FSJgenome'
PREFIX='AphCoe_V3_internal_Mar2024'

# Set working directory
cd ${TOPDIR}

# Run QUAST
quast.py -t 4 --eukaryote --large -o ${TOPDIR}/QUAST_${PREFIX} ${DRAFT}