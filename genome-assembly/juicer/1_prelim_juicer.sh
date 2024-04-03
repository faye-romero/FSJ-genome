#!/bin/bash
#SBATCH -p standard
#SBATCH --job-name="HifiSal_Dec22_prelim_juicer"
#SBATCH -N 1
#SBATCH -c 12
#SBATCH --mem=150G
#SBATCH --time=3-00:00:00

## Juicer, Part 1: Preliminary set-up ##
## Faye Romero. University of Rochester. 15 June 2023 ##

## Help pages:
# How to run Juicer with SLURM: https(:)//github(.)com/aidenlab/juicer/wiki/Running-Juicer-on-a-cluster
# Dependencies needed to run Juicer: https(:)//github(.)com/aidenlab/juicer/wiki/Installation#dependencies

# Make sure you've cloned the Juicer SLURM scripts directory from https(:)//github(.)com/aidenlab/juicer/tree/main/SLURM/scripts

## Load modules according to dependencies
module purge
module load circ
module load java
module load samtools/1.16.1
module load python/2.7.12
module load bwa/0.7.17
module load gnutls/3.5.15

## Set paths according to dependencies
export PATH=$HOME/bin:$PATH
export PATH=$PATH:/scratch/nchen11_lab/Juicer/scripts #path to Juicer Tools jar
export PATH=$PATH:/software/gnutls/3.5.15 #path to gnutls
export PATH=$PATH:/software/bwa/0.7.17 #path to bwa
export PATH=$PATH:/software/python/2.7.12 #path to python

## Define variables
DRAFT=/scratch/nchen11_lab/Juicer/references/HifiSal_Dec22.fasta # genome
PREFIX='HifiSal_Dec22' # genome name
ENZPOS_SCRIPT=/scratch/nchen11_lab/Juicer/scripts/generate_site_positions.py
JUICER_SCRIPT=/scratch/nchen11_lab/Juicer/scripts/juicer_Bluehive.sh # note that I've modified the main Juicer script for the UR cluster specifically. Please email me if you have any questions about my modifications.
RAW_HIC=/scratch/nchen11_lab/Juicer/fastq
JUICER_TOP_DIR=/scratch/nchen11_lab/Juicer
JUICER_SCRIPTS_DIR=/scratch/nchen11_lab/Juicer

# Set up directory structure
mkdir -p /scratch/nchen11_lab/Juicer
cd /scratch/nchen11_lab/Juicer
mkdir references restriction_sites fastq
ln -s $DRAFT ./references #soft-link genome into references folder

## Steps taken from p. 15 of Juicer assembly cookbook ##

#STEP 1. Index draft
cd ./references
bwa index $DRAFT
echo "bwa indexing done!"
cd ../

#STEP 2. Generate restriction site enzyme positions
python $ENZPOS_SCRIPT DpnII ./restriction_sites/${PREFIX} ./references/${PREFIX}.fasta
cd ./restriction_sites
echo "generate restriction sites file done!" date

#STEP 3. Generate chrom.sizes file
awk 'BEGIN{OFS="\t"}{print $1, $NF}' ${PREFIX}_DpnII.txt > ${PREFIX}.chrom.sizes
cd ../
echo "generate chrom.sizes file done!" date

#STEP 4. Run the juicer.sh alignment script using the 'juicer_terminal_command.txt' directly on the command line.
