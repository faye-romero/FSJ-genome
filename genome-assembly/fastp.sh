#!/bin/bash
#SBATCH -p rosalind
#SBATCH --job-name="fastp"
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=50G
#SBATCH --time=7-00:00:00
#SBATCH --output=fastp.03182023.log

## Fastp ##
## Faye Romero. University of Rochester. 18 March 2023 ##

# Load module
module load fastp/b1

# Make/set working directories
cd /scratch/nchen11_lab/Faye/ReferenceMarch23
mkdir ./fastp_out

# Set variables (paired-end reads for mom and dad)
DAD_ONE=/scratch/nchen11_lab/FSJ_WGS/rawData/H1172/H1172_CKDN220066516-1A_HT2VNDSX5_L1_1.fq.gz
DAD_TWO=/scratch/nchen11_lab/FSJ_WGS/rawData/H1172/H1172_CKDN220066516-1A_HT2VNDSX5_L1_2.fq.gz
MOM_ONE=/scratch/nchen11_lab/FSJ_WGS/rawData/H914/H914_CKDN230002560-1A_HT5NVDSX5_L1_1.fq.gz
MOM_TWO=/scratch/nchen11_lab/FSJ_WGS/rawData/H914/H914_CKDN230002560-1A_HT5NVDSX5_L1_2.fq.gz

# Run fastp
#usage: fastp -i in.R1.fq.gz -I in.R2.fq.gz -o out.R1.fq.gz -O out.R2.fq.gz
fastp -i $DAD_ONE -I $DAD_TWO -o ./fastp_out/fastp.H1172.R1.out.fq.gz -O ./fastp_out/fastp.H1172.R2.out.fq.gz
fastp -i $MOM_ONE -I $MOM_TWO -o ./fastp_out/fastp.H914.R1.out.fq.gz -O ./fastp_out/fastp.H914.R2.out.fq.gz
