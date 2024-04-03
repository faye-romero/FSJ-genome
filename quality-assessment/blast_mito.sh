#!/bin/bash
#SBATCH -p rosalind
#SBATCH --job-name="mito.blast"
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=100G
#SBATCH --time=5-00:00:00
#SBATCH --output=mito.blast.02142024.log

## BLAST for mitochondrial contamination ##
## Faye Romero. University of Rochester. 14 February 2023 ##

# Load modules and set PATH
module load blast
export PATH=/software/blast/2.10.0+/blastdb:$PATH

# Set wd and variables
cd /scratch/nchen11_lab/new_FSJgenome/mitochondria
QUERY='/scratch/nchen11_lab/new_FSJgenome/allmaps_out/HifiSal_Dec22_v2_LOD5update_MAP.fasta'
PREFIX='HifiSal_Dec22_v2_LOD5update_MAP'

# Florida Scrub-Jay mitochondria: NCBI accession NC_051467.1
echo "$(date).....blasting for the FSJ mitochondria......"
echo "$(date).....makeblastdb -in NC_051467.1.fasta -dbtype nucl -parse_seqids -out NC_051467.1.db....."
makeblastdb -in NC_051467.1.fasta -dbtype nucl -parse_seqids -out NC_051467.1.db

echo "$(date).....blastn -db NC_051467.1.db -query $GENOME -outfmt 6 -evalue 0.00001 -num_threads 8 -out ${PREFIX}.mito.blast.out.txt....."
blastn -db NC_051467.1.db -query $QUERY -outfmt 6 -evalue 0.00001 -num_threads 8 -out ${PREFIX}.mito.blast.out.txt
