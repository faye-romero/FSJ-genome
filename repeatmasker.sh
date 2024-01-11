#!/bin/bash
#SBATCH -p rosalind
#SBATCH --job-name="repeatmasker"
#SBATCH -N 1
#SBATCH -c 24
#SBATCH --mem=100G
#SBATCH --time=30-00:00:00
#SBATCH --output=repeatmasker.11302023.log
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=fromero3@ur.rochester.edu

## Script to run RepeatMasker on a genome assembly (with a local installation of RepeatMasker) ##
## Faye Romero, University of Rochester ##
## 30 Nov 2023 ##

# Set wd
cd /scratch/nchen11_lab/FSJgenome

# Set variables
LIBRARY='/scratch/nchen11_lab/FSJgenome/repeat_libraries/FSJv3_master_repeats.lib' # Repeat library. In this case, this is a combined repeat library which includes: custom repeat library from RepeatMasker, avian repeat library (from Valentina Peona), and the Steller's Jay curated repeat library "bCyaSte1.0.p_RepeatLibrary.fa"
GENOME='/scratch/nchen11_lab/FSJgenome/HifiSal_Dec22_v2_ALLMAPS.postyag_newnames.fasta'
HOME='/scratch/nchen11_lab/Repeat' # path to directory with RepeatModeler + RepeatMasker program folders
export PATH=/scratch/nchen11_lab/Repeat:$PATH # add HOME dir to your PATH

# Run RepeatMasker
echo "beginning. $(date)."
echo "${HOME}/RepeatMasker/RepeatMasker $GENOME \
-lib $LIBRARY \
-gff \
-s \
-a \
-xsmall \
-pa 24 \
-dir repeatmasker_out"

${HOME}/RepeatMasker/RepeatMasker $GENOME \
-lib $LIBRARY \
-gff \
-s \
-a \
-xsmall \
-pa 24 \
-dir repeatmasker_out

echo "finished. $(date)."