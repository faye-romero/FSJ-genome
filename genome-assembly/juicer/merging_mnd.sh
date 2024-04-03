#!/bin/bash
#SBATCH -p rosalind
#SBATCH --job-name="sorting_mnd"
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=50G
#SBATCH --time=3-00:00:00
#SBATCH --output=sorting_mnd.06142023.log
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=fromero3@ur.rochester.edu

## Combining replicate HiC libraries by combining merged_nodups.txt files

cd /scratch/nchen11_lab/Juicer

mnd1=/scratch/nchen11_lab/Juicer/aligned/merged_nodups.txt
mnd2=/scratch/nchen11_lab/Juicer_AltLib/aligned/merged_nodups.txt

#sort mnd files
#sort -k2,2d -k6,6d ${mnd1} > /scratch/nchen11_lab/Juicer/aligned/merged_nodups_sorted_402.txt
#sort -k2,2d -k6,6d ${mnd2} > /scratch/nchen11_lab/Juicer_AltLib/aligned/merged_nodups_sorted_403.txt

SORTED1=/scratch/nchen11_lab/Juicer/aligned/merged_nodups_sorted_402.txt
SORTED2=/scratch/nchen11_lab/Juicer_AltLib/aligned/merged_nodups_sorted_403.txt

#merge and sort mnd files
#sort -m ${SORTED1} ${SORTED2} > merged_nodups_sorted_comb.txt
sort -k2,2d -k6,6d merged_nodups_sorted_comb.txt > merged_nodups_sorted_comb_v2.txt

#next, generate a .hic file from the combined mnd file. see 'juicertoolspre.sh'
