#!/bin/bash
#SBATCH -p rosalind
#SBATCH --job-name="quast"
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=10G
#SBATCH --time=12:00:00
#SBATCH --output=quast.12282023.log
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=fromero3@ur.rochester.edu

## Command to run QUAST on a large eukaryotic genome ##
## Faye Romero, University of Rochester ##
## 11 Jan 2024 ##

# Load modules
module load python3/3.7.1
module load quast
module load gcc/9.1.0
module load gnuplot/5.2.0b2
module load gnutls/3.5.15
module load zlib/1.2.11

# Set variables
DRAFT='/scratch/nchen11_lab/FSJgenome/FSJv3_final_unloc.scaffs.fasta'
TOPDIR='/scratch/nchen11_lab/FSJgenome'
PREFIX='FSJv3_final_unloc.scaffs'

# Set wd
cd ${TOPDIR}

# Run QUAST
quast.py -t 4 --eukaryote --large -o ${TOPDIR}/quast_${PREFIX} ${DRAFT}