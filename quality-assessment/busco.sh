#!/bin/bash
#SBATCH -p rosalind
#SBATCH --job-name="busco"
#SBATCH -N 1
#SBATCH -c 24
#SBATCH --mem=150G
#SBATCH --time=7-00:00:00
#SBATCH --output=busco.03022024.log

## BUSCO ##
## Faye Romero. University of Rochester. 02 March 2024 ##

# Set working directory
wd='/scratch/nchen11_lab/new_FSJgenome'
cd $wd

# Load modules
module load busco/5.2.2

# Set variables
GENOME='/scratch/nchen11_lab/new_FSJgenome/FSJv3_internal_Feb2024_v2.fasta'

# Create and navigate into output directory
mkdir ${wd}/BUSCO-euk_FSJv3_internal_Feb2024_v2
cd ${wd}/BUSCO-euk_FSJv3_internal_Feb2024_v2

#Run BUSCO for eukaryota database
echo "$(date).....busco -i $GENOME --lineage_dataset eukaryota_odb10 -m genome -o BUSCO-eukaryota_FSJv3_internal_Feb2024_v2 -c 24"
busco -i $GENOME --lineage_dataset eukaryota_odb10 -m genome -o BUSCO-eukaryota_FSJv3_internal_Feb2024_v2 -c 24
cd $wd

# Create and navigate into output directory
mkdir ${wd}/BUSCO-aves_FSJv3_internal_Feb2024_v2
cd ${wd}/BUSCO-aves_FSJv3_internal_Feb2024_v2

#Run BUSCO for Aves database
echo "$(date).....busco -i $GENOME --lineage_dataset aves_odb10 -m genome -o BUSCO-aves_FSJv3_internal_Feb2024_v2 -c 24"
busco -i $GENOME --lineage_dataset aves_odb10 -m genome -o BUSCO-aves_FSJv3_internal_Feb2024_v2 -c 24
cd $wd