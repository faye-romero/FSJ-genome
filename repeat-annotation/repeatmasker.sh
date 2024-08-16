#!/bin/bash
#SBATCH -p rosalind
#SBATCH --job-name="repeatmasker"
#SBATCH -N 1
#SBATCH -c 21
#SBATCH --mem=100G
#SBATCH --nodelist=bhd0056
#SBATCH --time=7-00:00:00
#SBATCH --output=repeatmasker.07222024.log
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=fromero3@ur.rochester.edu

cd /scratch/nchen11_lab/Faye/please

GENOME='/scratch/nchen11_lab/Faye/please/FSJgenome_July2024_FINAL.fasta'
HOME='/scratch/nchen11_lab/Repeat'
export PATH=/scratch/nchen11_lab/Repeat:$PATH

#Set variables
LIBRARY='/scratch/nchen11_lab/Faye/please/repeatmodeler/FSJgenome_July2024_FINAL_master_repeats.lib'

#Run RepeatMasker
echo "$(date).....starting."
echo "${HOME}/RepeatMasker/RepeatMasker $GENOME \
-lib $LIBRARY \
-gff \
-s \
-a \
-xsmall \
-pa 21 \
-dir repeatmasker_out"

${HOME}/RepeatMasker/RepeatMasker $GENOME \
-lib $LIBRARY \
-gff \
-s \
-a \
-xsmall \
-pa 21 \
-dir repeatmasker_out

echo "$(date).....finished."

#-gff = outputs a .gff file
#-s = slow search
#-a = shows the alignments in a .align output file
#-xsmall = returns repetitive regions in lowercase (rest capitals) rather than masked
#-pa = number of processors to use in parallel

#Output = RepeatMasker returns a .masked file containing the query sequence(s) with all identified repeats and low complexity sequences masked. These masked sequences are listed and annotated in the .out file. The masked sequences are printed in the same order as they are in the submitted file, whereas the sequences are presented alphabetically in the annotation table. The .tbl file is a summary of the repeat content of the analyzed sequence.
