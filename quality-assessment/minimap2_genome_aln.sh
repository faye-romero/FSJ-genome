#!/bin/bash
#SBATCH -p rosalind
#SBATCH --job-name="minimap2.genome.aln"
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=50G
#SBATCH --time=3-00:00:00
#SBATCH --output=minimap2.genome.aln.log

## Minimap2 interspecies whole genome alignment ##
## Faye Romero. University of Rochester. 15 March 2024 ##

# Load modules
module load minimap2/2.26
module load samtools/1.9

# Set variables
REF='/scratch/nchen11_lab/Faye/AvesReferenceGenomes/bCorMon1.pri_genomic_clean.fna'
RPREFIX='bCorMon1.pri_genomic_clean'
QUERY='/scratch/nchen11_lab/new_FSJgenome/allmaps_out_v2/FSJv3_internal_Feb2024.fasta'
QPREFIX='FSJv3_internal_Feb2024'

# Set wd
wd='/scratch/nchen11_lab/new_FSJgenome/whole_genome_alns'
cd $wd

# Run minimap2. -x asm20 indicates relatively diverged species.
echo "$(date).....minimap2 -x asm20 -t 8 $REF $QUERY > ${RPREFIX}-VERSUS-${QPREFIX}.paf"
minimap2 -x asm20 -t 8 $REF $QUERY > ${RPREFIX}-VERSUS-${QPREFIX}.paf