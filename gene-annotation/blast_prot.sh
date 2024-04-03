#!/bin/bash
#SBATCH -p rosalind
#SBATCH --job-name="prot.blast"
#SBATCH -N 1
#SBATCH -c 24
#SBATCH --mem=350G
#SBATCH --time=7-00:00:00
#SBATCH --output=prot.blast.03272024.log

## Protein BLAST ##
## Faye Romero. University of Rochester. 27 March 2024 ##

# Load modules and set PATH
module load blast
export PATH=/software/blast/2.10.0+/blastdb:$PATH

# Set variables and softlinks
MYPROTEINS='/scratch/nchen11_lab/FSJ_GeneAnnotation/braker.clean.March2024.aa'
NAME='braker.clean.March2024.aa'
ln -s $MYPROTEINS /scratch/nchen11_lab/uniprot_db

cd /scratch/nchen11_lab/uniprot_db

# Blast against Swiss-Prot protein database
# https(:)//www.uniprot(.)org/uniprotkb?query=reviewed:true

echo "$(date).....creating blast db from Swiss-Prot fasta sequence....."
makeblastdb -in uniprot_sprot_11Dec2023.fasta -dbtype prot -out uniprot_sprot_11Dec2023.db

echo "$(date).....blasting BRAKER3 predicted protein sequences against Swiss-Prot protein database....."
echo "blastp -db uniprot_sprot_11Dec2023.db -query $MYPROTEINS -evalue 0.000001 -outfmt 6 -num_threads 24 -out braker.aa_uniprot_sprot_27March2024.out.blast"
blastp -db uniprot_sprot_11Dec2023.db -query $NAME -evalue 0.000001 -outfmt 6 -num_threads 24 -out braker.aa_uniprot_sprot_27March2024.out.blast