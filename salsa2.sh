#!/bin/bash
#SBATCH -p standard
#SBATCH --job-name="salsa"
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=300G
#SBATCH --time=3-00:00:00
#SBATCH --output=salsa.09112023.log
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=fromero3@ur.rochester.edu

## Script to run SALSA2 scaffolder ##
## Faye Romero, University of Rochester ##
## 11 Sept 2023 ##

# Set wd
cd /scratch/nchen11_lab/Faye/ReferenceAug23

# Load modules
module load bedtools/2.30.0
module load samtools/1.9
module load salsa/2.3

# Set variables
IN='HiC_to_clean_Hifiasm_dip.hap1.H1172_MAP' # prefix of input contig file
IN_DIR='/scratch/nchen11_lab/Faye/ReferenceAug23/arima_out/arima_out_H1172post' # directory of input contig file
CONTIGS='/scratch/nchen11_lab/Faye/ReferenceAug23/allmaps_out_08302023/hap1.H1172_out/clean_Hifiasm_dip.hap1.H1172_MAP.fasta' # input contig file
OUT_DIR='/scratch/nchen11_lab/Faye/ReferenceAug23/salsa_out' # directory of output scaffold file
OUT='Hifiasm_dip.hap1.H1172_MAP_salsa' # prefix of output scaffold file
SALSAPIPE='/software/salsa/2.3/run_pipeline.py' # path to salsa pipeline

# Convert bam to bed format
bamToBed -i ${IN_DIR}/${IN}.bam > ${IN_DIR}/${IN}.bed
sort -k 4 ${IN_DIR}/${IN}.bed > tmp && mv tmp ${IN_DIR}/${IN}.bed

# Generate contig lengths file with samtools faidx 
samtools faidx $CONTIGS

# Run SALSA2
python $SALSAPIPE -a $CONTIGS -l ${CONTIGS}.fai -b ${IN_DIR}/${IN}.bed -e GATC -o ${OUT_DIR}/${OUT} -i 5 -m yes
