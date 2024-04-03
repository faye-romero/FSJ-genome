#!/bin/bash
#SBATCH -p rosalind
#SBATCH --job-name="Arima_merge"
#SBATCH -N 1
#SBATCH -c 24
#SBATCH --mem=320G
#SBATCH --time=3-00:00:00

## Arima Mapping Pipeline, Part 3 (to merge Arima output from replicate Hi-C libraries)##
## Faye Romero. University of Rochester. 23 March 2023 ##

# Set wd
cd /scratch/nchen11_lab/Faye/genomes/NewReference

# Load modules
module purge
module load perl
module load bwa/0.7.17
module load samtools/1.9
module load picard/2.12.0

##############################################
# ARIMA GENOMICS MAPPING PIPELINE 02/08/2019 #
##############################################

#Below find the commands used to map HiC data.

#Replace the variables at the top with the correct paths for the locations of files/programs on your system.

#This bash script will map one paired end HiC dataset (read1 & read2 fastqs). Feel to modify and multiplex as you see fit to work with your volume of samples and system.

##########################################
# Commands #
##########################################

LABEL='HiCtonewHifiasm_revcom'
IN_DIR='/scratch/nchen11_lab/hic/hic_raw'
REF='/scratch/nchen11_lab/Faye/genomes/NewReference/newHifiasm_revcom.fa'
OUT_DIR='/scratch/nchen11_lab/Faye/genomes/NewReference/HiC'
FAIDX='/scratch/nchen11_lab/Faye/genomes/NewReference/newHifiasm_revcom.fa.fai'
FILTER='/scratch/nchen11_lab/Arima_scripts/filter_five_end.pl'
COMBINER='/scratch/nchen11_lab/Arima_scripts/two_read_bam_combiner.pl'
STATS='/scratch/nchen11_lab/Arima_scripts/get_stats.pl'
PICARD='/software/picard/2.12.0/picard.jar'
TMP_DIR='/scratch/nchen11_lab/Faye/genomes/NewReference/HiC/tmp'
REP_LABEL=$LABEL\_rep2
MAPQ_FILTER=10
CPU=24

#########################################################################################################################################
###                                       How to Accommodate Biological Replicates                                                    ###
### This pipeline is currently built for processing a single sample with one read1 and read2 fastq file.                              ###
### Biological replicates (eg. multiple libraries made from the same sample) should be merged before proceeding with subsequent steps.###
### The code below is an example of how to merge biological replicates.                                                               ###
#########################################################################################################################################
echo "### Step 5: Merge libraries"

INPUTS_BIOLOGICAL_REPS=('INPUT=/scratch/nchen11_lab/Faye/genomes/NewReference/HiC/HiCtonewHifiasm_revcom_rep1.bam' 'INPUT=/scratch/nchen11_lab/Faye/genomes/NewReference/HiC/HiCtonewHifiasm_revcom_rep2.bam') #BAM files you want combined as biological replicates

java -Xmx8G -Djava.io.tmpdir=temp/ -jar $PICARD MergeSamFiles $INPUTS_BIOLOGICAL_REPS OUTPUT=$OUT_DIR/$LABEL.bam USE_THREADING=TRUE ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT

samtools index $OUT_DIR/$LABEL.bam

perl $STATS $OUT_DIR/$LABEL.bam > $OUT_DIR/$LABEL.bam.stats

echo "Finished Mapping Pipeline through merging Biological Replicates"
