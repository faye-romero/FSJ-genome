#!/bin/bash
#SBATCH -p rosalind
#SBATCH --job-name="nt.blast"
#SBATCH -N 1
#SBATCH -c 32
#SBATCH --mem=350G
#SBATCH --time=5-00:00:00
#SBATCH --output=nt.blast.02162024.log

## BLAST against nt database (for Blobtoolkit decontamination screening) ##
## Faye Romero. University of Rochester. 20 February 2024 ##

# Load modules
module load blast

# Set wd
cd /scratch/nchen11_lab/new_FSJgenome

# Set variables and BLASTDB
QUERY='/scratch/nchen11_lab/new_FSJgenome/allmaps_out/HifiSal_Dec22_v2_LOD5update_MAP.fasta'
PREFIX='HifiSal_Dec22_v2_LOD5update_MAP'
export BLASTDB=/software/blast/2.10.0+/blastdb:/scratch/fromero3/taxdb # add appropriate databases to BLASTDB

echo "$(date).....blasting against nt database....."
echo "blastn -db /software/blast/2.10.0+/blastdb/nt \
-query $QUERY \
-outfmt 6 qseqid staxids sskingdoms sscinames scomnames bitscore std \
-max_target_seqs 10 \
-max_hsps 1 \
-evalue 1e-25 \
-num_threads 32 \
-out ${PREFIX}.nt.blast.out.txt"

blastn -db /software/blast/2.10.0+/blastdb/nt \
-query $QUERY \
-outfmt "6 qseqid staxids sskingdoms sscinames scomnames bitscore std" \
-max_target_seqs 10 \
-max_hsps 1 \
-evalue 1e-25 \
-num_threads 32 \
-out ${PREFIX}.nt.blast.out.txt

echo "$(date).....done."