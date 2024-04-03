#!/bin/bash
#SBATCH -p rosalind
#SBATCH --job-name="blobtools"
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=100G
#SBATCH --time=7-00:00:00
#SBATCH --output=blobtools.02202024.log

## Blobtoolkit (decontamination screening) ##
## Faye Romero. University of Rochester. 20 February 2024 ##

# Load conda env
module load miniconda3
source activate blobtoolkit

# Set wd
wd='/scratch/nchen11_lab/new_FSJgenome'
cd $wd

# Set variables
GENOME='/scratch/nchen11_lab/new_FSJgenome/allmaps_out/HifiSal_Dec22_v2_LOD5update_MAP.fasta'
PREFIX='HifiSal_Dec22_v2_LOD5update_MAP'

## Create blob object ##
echo '$(date)... blobtools create --fasta $GENOME FSJblob'
blobtools create --fasta $GENOME FSJblob

## Add BUSCO scores ##
#To get BUSCO results, see busco.sh

echo '$(date).....blobtools add --busco ${wd}/BUSCO-euk_HifiSal_Dec22_v2_LOD5update_MAP/BUSCO-eukaryota_HifiSal_Dec22_v2_LOD5update_MAP/run_eukaryota_odb10/full_table.tsv'
blobtools add --busco ${wd}/BUSCO-euk_HifiSal_Dec22_v2_LOD5update_MAP/BUSCO-eukaryota_HifiSal_Dec22_v2_LOD5update_MAP/run_eukaryota_odb10/full_table.tsv

echo '$(date).....blobtools add --busco ${wd}/BUSCO-aves_HifiSal_Dec22_v2_LOD5update_MAP/BUSCO-aves_HifiSal_Dec22_v2_LOD5update_MAP/run_aves_odb10/full_table.tsv'
blobtools add --busco ${wd}/BUSCO-aves_HifiSal_Dec22_v2_LOD5update_MAP/BUSCO-aves_HifiSal_Dec22_v2_LOD5update_MAP/run_aves_odb10/full_table.tsv

## Add BLAST hits ##
# To get BLAST results, see blast_nt.sh
# BLAST hits: ./HifiSal_Dec22_v2_ALLMAPS.postyag_newnames.blast.out

echo '$(date)... blobtools add --hits HifiSal_Dec22_v2_LOD5update_MAP.nt.blast.out.txt --taxrule bestsumorder --taxdump /scratch/fromero3/taxdump --threads 8 FSJblob'
blobtools add --hits HifiSal_Dec22_v2_LOD5update_MAP.nt.blast.out \
--taxrule bestsumorder \
--taxdump /scratch/fromero3/taxdump \
--hits-cols 1=qseqid,2=staxids,6=bitscore,8=sseqid,13=qstart,14=qend,17=evalue \
FSJblob

## Add read coverage ##
#To get read coverage, see minimap2_btk.sh

echo '$(date)... blobtools add --cov ${wd}/HifiSal_Dec22_v2_LOD5update_MAP.reads.bam --threads 8 FSJblob' 
blobtools add --cov ${wd}/HifiSal_Dec22_v2_LOD5update_MAP.reads.bam --threads 8 FSJblob