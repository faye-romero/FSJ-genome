#!/bin/bash
#SBATCH -p rosalind
#SBATCH --job-name="repeatmodeler.v2"
#SBATCH -N 1
#SBATCH -c 21
#SBATCH --mem=150G
#SBATCH --time=7-00:00:00
#SBATCH --output=repeatmodeler.03082024.log

## RepeatModeler ##
## Faye Romero. University of Rochester. 08 March 2024 ##

# Load modules
module unload python python3 miniconda miniconda3
module load miniconda3
source activate tesorter # Activate conda env with tesorter installed
module load perl

# Set wd
wd='/scratch/nchen11_lab/new_FSJgenome/repeatmodeler'
cd $wd

# Set variables
GENOME='/scratch/nchen11_lab/new_FSJgenome/AphCoe_V3_internal_Mar2024.fasta'
HOME='/scratch/nchen11_lab/Repeat'

# Set paths to software. I am using a local installation of RepeatModeler/Masker.
export PATH=/scratch/nchen11_lab/Repeat:$PATH
export PATH=/scratch/nchen11_lab/hmmer-3.4/bin:$PATH
export PATH=/scratch/nchen11_lab/Repeat/RepeatMasker:$PATH
export PATH=/scratch/nchen11_lab/Repeat/RECON-1.08/bin:$PATH
export PATH=/scratch/nchen11_lab/Repeat/RepeatScout-1.0.6:$PATH
export PATH=/scratch/nchen11_lab/Repeat/cd-hit-v4.8.1-2019-0228:$PATH
export PATH=/scratch/nchen11_lab/Repeat/RepeatModeler-2.0.4:$PATH
export PATH=/scratch/nchen11_lab/Repeat/rmblast-2.13.0/bin:$PATH
export PATH=/scratch/nchen11_lab/Repeat/genometools-1.6.2/bin:$PATH
export PATH=/scratch/nchen11_lab/Repeat/LTR_retriever-2.9.0:$PATH
export PATH=/scratch/nchen11_lab/Repeat/mafft-7.490-with-extensions/bin:$PATH
export PATH=/scratch/nchen11_lab/Repeat/NINJA-0.95-cluster_only/NINJA:$PATH

# Build the database
DB='AphCoe_V3_internal_Mar2024_db' #Set database name

echo "building RM database. $(date)"
echo "${HOME}/RepeatModeler-2.0.4/BuildDatabase -name $DB $GENOME"
${HOME}/RepeatModeler-2.0.4/BuildDatabase -name $DB $GENOME

echo "finished building RM database. $(date)"

#Run RepeatModeler
cd $wd

echo "running RepeatModeler within the tesorter conda env. $(date)"
echo "${HOME}/RepeatModeler-2.0.4/RepeatModeler -database $DB \
-engine ncbi \
-threads 21 \
-LTRStruct"

${HOME}/RepeatModeler-2.0.4/RepeatModeler -database $DB \
-threads 21 \
-LTRStruct

echo "finished running RepeatModeler. $(date)"
