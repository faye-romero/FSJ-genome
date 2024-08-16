#!/bin/bash
#SBATCH -p rosalind
#SBATCH --job-name="cdhitest"
#SBATCH -N 1
#SBATCH -c 2
#SBATCH --mem=50G
#SBATCH --time=3-00:00:00
#SBATCH --output=cdhitest.log

## Script to run cd-hit-test on RepeatModeler library output ##
## Faye Romero. University of Rochester. ##
## 18 April 2024 ##

# Set working directory of your choosing
wd='/scratch/nchen11_lab/Faye/please/repeatmodeler'
cd $wd

# Set variables
HOME='/scratch/nchen11_lab/Repeat' # location of local RepeatModeler/Masker installs
RMLIB='FSJgenome_July2024_FINAL_db-families_oneline'

# Run cd-hit-est
echo "$(date)....running cd-hit-est to reduce the redundancy of the RepeatModeler output by clustering sequences that have more than 80% sequence identity....."
echo ".....${HOME}/cd-hit-est -i ${RMLIB}.fa -o ${RMLIB}.reduced.fa -d 0 -aS 0.8 -c 0.8 -G 0 -g 1 -b 500....."

${HOME}/cd-hit-v4.8.1-2019-0228/cd-hit-est -i ${RMLIB}.fa -o ${RMLIB}.reduced.fa -d 0 -aS 0.8 -c 0.8 -G 0 -g 1 -b 500

echo "$(date).....finished....."

# Parameters:
# “-d” specifies the length of the fasta header (with 0 meaning “keep everything up until the first space)
# “-aS” defines the alignment coverage for the shorter sequence, when set to 0.8 the shorter sequence must have at least 80% of its length matched to the representative sequence, to comply with the 80-80-80 rule
# -c sets the global sequence identity, which is set to 80% to comply with the 80-80-80 rule
# -G defines how identity is calculated, set it to “0” for local sequence alignment and “1” for global alignment. Local (G = 0) ensures that the 80% identity threshold applies only to the aligned stretches between the queried and representative sequences instead of the complete end-to-end (global) alignment of the two.
# -g chooses the clustering algorithm: each sequence is clustered to either first (0) or best (1) representative sequence. 
# -b defines the alignment’s bandwidth and has been set to 500. This setting avoids splitting sequences into different clusters where an insertion/deletion causes a break in the synteny of matching identities which, without the indel, is at least 80% of the sequence. 