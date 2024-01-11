#!/bin/bash
#SBATCH -p rosalind
#SBATCH --job-name="repeatmodeler"
#SBATCH -N 1
#SBATCH -c 21
#SBATCH --mem=100G
#SBATCH --time=30-00:00:00
#SBATCH --output=repeatmodeler.11132023.log
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=fromero3@ur.rochester.edu

## Script to run RepeatModeler on a genome assembly (with a local installation of RepeatModeler) ##
## Faye Romero ##
## 13 Nov 2023 ##

# Set wd
cd /scratch/nchen11_lab/FSJgenome

# Set variables
GENOME='/scratch/nchen11_lab/FSJgenome/HifiSal_Dec22_v2_ALLMAPS.postyag_newnames.fasta'
HOME='/scratch/nchen11_lab/Repeat' # path to directory with RepeatModeler + RepeatMasker program folders
export PATH=/scratch/nchen11_lab/Repeat:$PATH # add HOME dir to your PATH

# Build the database
echo "building RM database. $(date)"
echo "${HOME}/RepeatModeler-2.0.4/BuildDatabase -name FSJv3_db -engine ncbi $GENOME"

${HOME}/RepeatModeler-2.0.4/BuildDatabase -name FSJv3_db $GENOME

echo "finished building RM database. $(date)"

# Run RepeatModeler
echo "running RepeatModeler. $(date)"
echo "${HOME}/RepeatModeler-2.0.4/RepeatModeler -database FSJv3_db \
-engine ncbi \
-threads 21 \
-LTRStruct"

${HOME}/RepeatModeler-2.0.4/RepeatModeler -database FSJv3_db \
-threads 21 \
-LTRStruct

echo "finished running RepeatModeler. $(date)"
