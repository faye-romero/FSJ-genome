#!/bin/bash
#SBATCH -p rosalind
#SBATCH --job-name="repeatmasker.FSJV3"
#SBATCH -N 1
#SBATCH -c 21
#SBATCH --mem=100G
#SBATCH --time=7-00:00:00
#SBATCH --output=repeatmasker.V3.03152024.log

## RepeatMasker ##
## Faye Romero. University of Rochester. 08 March 2024 ##

# Set wd
cd /scratch/nchen11_lab/new_FSJgenome

# Set variables and PATHs. I am using a local installation of RepeatModeler/Masker.
GENOME='/scratch/nchen11_lab/new_FSJgenome/AphCoe_V3_internal_Mar2024.fasta'
HOME='/scratch/nchen11_lab/Repeat'
export PATH=/scratch/nchen11_lab/Repeat:$PATH

#Set variables
LIBRARY='/scratch/nchen11_lab/new_FSJgenome/repeatmodeler/AphCoe_V3_internal_Mar2024_master_repeatlib.fa' # This is a concatenated lib of de novo Florida Scrub-Jay repeats from RepeatModeler, Valentina Peona's 2021 bird repeat library, and a Steller's Jay repeat library (Benham et al. 2023).

#Run RepeatMasker
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

#-gff = outputs a .gff file
#-s = slow search
#-a = shows the alignments in a .align output file
#-xsmall = returns repetitive regions in lowercase (rest capitals) rather than masked
#-pa = number of processors to use in parallel

#Output = RepeatMasker returns a .masked file containing the query sequence(s) with all identified repeats and low complexity sequences masked. These masked sequences are listed and annotated in the .out file. The masked sequences are printed in the same order as they are in the submitted file, whereas the sequences are presented alphabetically in the annotation table. The .tbl file is a summary of the repeat content of the analyzed sequence.
