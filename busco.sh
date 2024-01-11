#!/bin/bash
#SBATCH -p rosalind
#SBATCH --job-name="busco"
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=150G
#SBATCH --time=3-00:00:00
#SBATCH --output=busco.12282023.log
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=fromero3@ur.rochester.edu

## Command to run BUSCO on a genome assembly, with auto lineage select ##
## Faye Romero, University of Rochester ##
## 11 Jan 2024 ##

# Set wd
cd /scratch/nchen11_lab/FSJgenome

# Load modules
module load busco/5.2.2

# Set variables
GENOME='/scratch/nchen11_lab/FSJgenome/FSJv3_final_unloc.scaffs.fasta'

# Run BUSCO
busco -i $GENOME --auto-lineage-euk -m genome -o BUSCO-auto_FSJv3_final_unloc.scaffs -c 8

