#!/bin/bash
#SBATCH -p rosalind
#SBATCH --job-name="visualize_juicer"
#SBATCH -N 1
#SBATCH -c 2
#SBATCH --mem=300G
#SBATCH --time=3-00:00:00
#SBATCH --output=visualize_juicer.HifiSal_Dec22_comblib.06152023.log

## Juicer, Part 3: Generate HiC file ##
## Faye Romero. University of Rochester. 15 June 2023 ##

# Set wd
cd /scratch/nchen11_lab/Juicer

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

# Set variables
FASTA=/scratch/nchen11_lab/Juicer/references/HifiSal_Dec22.fasta # genome
MNDCOMB=/scratch/nchen11_lab/Juicer/merged_nodups_sorted_comb_v2.txt #merged no dupes file (combined from replicate libraries)
PREFIX="HifiSal_Dec22_comblib" # output name
JUICER_SCRIPTS_DIR=/scratch/nchen11_lab/Juicer

mkdir ./visualize

#Step 3. convert your reference fasta file to a .assembly file using generate-sorted-assembly-file-from-fasta.awk
awk -f ${JUICER_SCRIPTS_DIR}/scripts/3d-dna-master/utils/generate-sorted-assembly-file-from-fasta.awk $FASTA > ./visualize/$PREFIX.assembly

#Step 4. Use the merged_nodups.txt file and the .assembly file to run run-assembly-visualizer.sh.
#./run-assembly-visualizer.sh [options] <path_to_input_assembly_file> <path_to_input-mnd-file>.
cd ${JUICER_SCRIPTS_DIR}/visualize
bash ${JUICER_SCRIPTS_DIR}/scripts/3d-dna-master/visualize/run-assembly-visualizer-Bluehive.sh $PREFIX.assembly $MNDCOMB

# note that I've modified the run-assembly-visualizer.sh script for the UR cluster specifically. Please email me if you have any questions about my modifications.

# Now, you can load the .assembly file and the .hic file into Juicebox for visualization/manual curation. This is step is "Part 4".