#!/bin/bash
#SBATCH -p standard
#SBATCH --job-name="salsa"
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=300G
#SBATCH --time=3-00:00:00
#SBATCH --output=salsa.09112023.log

## SALSA ##
## Faye Romero. University of Rochester. 09 September 2023 ##

# Set wd
cd /scratch/nchen11_lab/Faye/ReferenceAug23

# Load modules
module load bedtools/2.30.0
module load samtools/1.9
module load salsa/2.3

# Set variables
IN='HiC_to_clean_Hifiasm_dip.hap1.H1172_MAP' # Prefix of your genome
IN_DIR='/scratch/nchen11_lab/Faye/ReferenceAug23/arima_out/arima_out_H1172post' # Dir of your finished Arima mapping
CONTIGS='/scratch/nchen11_lab/Faye/ReferenceAug23/allmaps_out_08302023/hap1.H1172_out/clean_Hifiasm_dip.hap1.H1172_MAP.fasta' # Genome
OUT_DIR='/scratch/nchen11_lab/Faye/ReferenceAug23/salsa_out' # Out dir
OUT='Hifiasm_dip.hap1.H1172_MAP_salsa' # Out prefix
SALSAPIPE='/software/salsa/2.3/run_pipeline.py' # Path to SALSA pipeline

# convert bam to bed format
bamToBed -i ${IN_DIR}/${IN}.bam > ${IN_DIR}/${IN}.bed
sort -k 4 ${IN_DIR}/${IN}.bed > tmp && mv tmp ${IN_DIR}/${IN}.bed

# generate contig lengths file
samtools faidx $CONTIGS

# run salsa2
python $SALSAPIPE -a $CONTIGS -l ${CONTIGS}.fai -b ${IN_DIR}/${IN}.bed -e GATC -o ${OUT_DIR}/${OUT} -i 5 -m yes
