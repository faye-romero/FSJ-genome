#!/bin/bash
#SBATCH -p rosalind
#SBATCH --job-name="minimap2.read.aln"
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=150G
#SBATCH --time=30-00:00:00
#SBATCH --output=depth.log

## Get genome-wide read coverage using minimap2 (for Blobtoolkit decontamination screening) ##
## Faye Romero. University of Rochester. 20 February 2024 ##

# Load modules
module load minimap2/2.26
module load samtools/1.9

# Set variables
GENOME='/scratch/nchen11_lab/new_FSJgenome/allmaps_out/HifiSal_Dec22_v2_LOD5update_MAP.fasta'
READS='/scratch/nchen11_lab/pacbio/hifi_reads_comb.fastq.gz' # raw PacBio reads
PREFIX='HifiSal_Dec22_v2_LOD5update_MAP'

# Set wd
wd='/scratch/nchen11_lab/new_FSJgenome'
cd $wd

echo "$(date).....mapping raw reads to genome HifiSal_Dec22_v2_LOD5update_MAP for btk."
echo "#SBATCH -N 1 -c 8 --mem=150G"
echo "minimap2 -ax map-hifi -t 8 $GENOME $READS | samtools sort -@ 8 -O bam -o HifiSal_Dec22_v2_LOD5update_MAP.reads.bam"

minimap2 -ax map-hifi -t 8 $GENOME $READS | samtools sort -@ 8 -O bam -o HifiSal_Dec22_v2_LOD5update_MAP.reads.bam

echo "$(date).....done."