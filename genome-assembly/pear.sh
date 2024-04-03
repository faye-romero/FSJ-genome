#!/bin/bash
#SBATCH -p rosalind
#SBATCH --job-name="pear"
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --mem=200G
#SBATCH --time=7-00:00:00
#SBATCH --output=pear.03182023.log

## PEAR ##
## Faye Romero. University of Rochester. 18 March 2023 ##

# Load module
module load pear/0.9.11

# Set/make working directories
cd /scratch/nchen11_lab/Faye/ReferenceMarch23/
mkdir ./pear_out

# Unzip fqs, if necessary
#gunzip /scratch/nchen11_lab/Faye/ReferenceMarch23/fastp_out/fastp.H1172.R1.out.fq.gz
#gunzip /scratch/nchen11_lab/Faye/ReferenceMarch23/fastp_out/fastp.H1172.R2.out.fq.gz
#gunzip /scratch/nchen11_lab/Faye/ReferenceMarch23/fastp_out/fastp.H914.R1.out.fq.gz
#gunzip /scratch/nchen11_lab/Faye/ReferenceMarch23/fastp_out/fastp.H914.R2.out.fq.gz

# Set variables (clean paired-end reads for mom and dad)
DAD_ONE=/scratch/nchen11_lab/Faye/ReferenceMarch23/fastp_out/fastp.H1172.R1.out.fq
DAD_TWO=/scratch/nchen11_lab/Faye/ReferenceMarch23/fastp_out/fastp.H1172.R2.out.fq
MOM_ONE=/scratch/nchen11_lab/Faye/ReferenceMarch23/fastp_out/fastp.H914.R1.out.fq
MOM_TWO=/scratch/nchen11_lab/Faye/ReferenceMarch23/fastp_out/fastp.H914.R2.out.fq

#usage: pear -f forward_read.fastq -r reverse_read.fastq -o output_prefix
pear -f $DAD_ONE -r $DAD_TWO -o ./pear_out/cleaned.merged.H1172 -j 16
pear -f $MOM_ONE -r $MOM_TWO -o ./pear_out/cleaned.merged.H914 -j 16

# Re-zip fqs
gzip /scratch/nchen11_lab/Faye/ReferenceMarch23/fastp_out/fastp.H1172.R1.out.fq
gzip /scratch/nchen11_lab/Faye/ReferenceMarch23/fastp_out/fastp.H1172.R2.out.fq
gzip /scratch/nchen11_lab/Faye/ReferenceMarch23/fastp_out/fastp.H914.R1.out.fq
gzip /scratch/nchen11_lab/Faye/ReferenceMarch23/fastp_out/fastp.H914.R2.out.fq
