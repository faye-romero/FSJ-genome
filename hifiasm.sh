#!/bin/bash
#SBATCH -p rosalind
#SBATCH --job-name="hifiasm"
#SBATCH -N 1
#SBATCH -c 32
#SBATCH --mem=350G
#SBATCH --time=3-00:00:00
#SBATCH --output=hifiasm.log
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=fromero3@ur.rochester.edu

## Script to run Hifiasm on raw reads ##
## Faye Romero, University of Rochester ##
## 6 May 2022 ##

# Set wd
cd /scratch/nchen11_lab/FSJgenome/HiFi

# Load modules
module load hifiasm/0.16.1

# Run Hifiasm
hifiasm -t 32 -o FSJ.hifiCut.3May2022.asm FSJ.hifi.cut.fastq.gz

# Convert outputs of Hifiasm into FASTA files
# awk '/^S/{print ">"$2;print $3}' FSJ.hifiCut.3May2022.asm.bp.p_ctg.gfa > FSJ3May22.p_ctg.fa
# awk '/^S/{print ">"$2;print $3}' FSJ.hifiCut.3May2022.asm.bp.hap1.p_ctg.gfa > FSJ3May22.hap1.p_ctg.fa
# awk '/^S/{print ">"$2;print $3}' FSJ.hifiCut.3May2022.asm.bp.hap2.p_ctg.gfa > FSJ3May22.hap2.p_ctg.fa