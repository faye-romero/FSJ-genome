#!/bin/bash
#SBATCH -p rosalind
#SBATCH --job-name="post_curation_juicer"
#SBATCH -N 1
#SBATCH -c 12
#SBATCH --mem=200G
#SBATCH --time=30-00:00:00
#SBATCH --output=post_curation_juicer.07202023.log
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=fromero3@ur.rochester.edu

## Juicer, Part 5: Generate a new genome from your Juicer manual curation ##
## Faye Romero. University of Rochester. 15 June 2023 ##

# Load modules
module purge
module load circ
module load jre/1.8.0_111
module load cuda/11.3
module load gcc/11.2.0
module load samtools/1.16.1
module load python/2.7.12
module load bwa/0.7.17
module load gnutls/3.5.15
module load parallel/20180222

# Set PATH variable
export PATH=$HOME/bin:$PATH
export PATH=$PATH:/scratch/nchen11_lab/Juicer/scripts
export PATH=$PATH:/software/gnutls/3.5.15
export PATH=$PATH:/software/gcc/11.2.0
export PATH=$PATH:/software/bwa/0.7.17
export PATH=$PATH:/software/python/2.7.12
export PATH=$PATH:/software/samtools/1.16.1
export PATH=$PATH:/software/cuda/11.3
export PATH=$PATH:/software/jre/1.8.0_111
export PATH=$PATH:/software/parallel/20180222

# Set wd
cd /scratch/nchen11_lab/Juicer

# Set variables
CURATED=/scratch/nchen11_lab/Juicer/postcuration_July2023/HifiSal_Dec22_comblib.review.07172023.assembly # Exported from JuiceBox, after manual curation is finished
FASTA=/scratch/nchen11_lab/Juicer/references/HifiSal_Dec22.fasta # genome
MND=/scratch/nchen11_lab/Juicer/merged_nodups_sorted_comb_v2.txt # merged_nodups.txt file
JUICER_TOP_DIR=/scratch/nchen11_lab/Juicer

# Set output dir
mkdir -p $JUICER_TOP_DIR/post_curation_July2023
cd $JUICER_TOP_DIR/postcuration_July2023

# Run script
sh /scratch/nchen11_lab/Juicer/scripts/3d-dna-master/run-asm-pipeline-post-review.sh --sort-output -r $CURATED $FASTA $MND