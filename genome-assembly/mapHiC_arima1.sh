#!/bin/bash
#SBATCH -p rosalind
#SBATCH --job-name="Arima_1"
#SBATCH -N 1
#SBATCH -c 24
#SBATCH --mem=320G
#SBATCH --time=4-00:00:00
#SBATCH --output=Arima1_03232023.log

## Arima Mapping Pipeline, Part 1 (for first HiC library)##
## Faye Romero. University of Rochester. 23 March 2023 ##

# Set working directory
cd /scratch/nchen11_lab/Faye/ReferenceMarch23

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

LABEL='HiCtoHifiasm_prim'
IN_DIR='/scratch/nchen11_lab/hic/hic_raw'
REF='/scratch/nchen11_lab/Faye/ReferenceMarch23/hifiasm_triobin_out_03202023/Hifiasm_prim.fa'
OUT_DIR='/scratch/nchen11_lab/Faye/ReferenceMarch23/arima_out'
FAIDX='/scratch/nchen11_lab/Faye/ReferenceMarch23/hifiasm_triobin_out_03202023/Hifiasm_prim.fa.fai'
FILTER='/scratch/nchen11_lab/Arima_scripts/filter_five_end.pl'
COMBINER='/scratch/nchen11_lab/Arima_scripts/two_read_bam_combiner.pl'
STATS='/scratch/nchen11_lab/Arima_scripts/get_stats.pl'
PICARD='/software/picard/2.12.0/picard.jar'
TMP_DIR='/scratch/nchen11_lab/Faye/ReferenceMarch23/arima_out/tmp'
REP_LABEL=$LABEL\_rep1
MAPQ_FILTER=10
CPU=24

echo "### Step 0: Check output directories exist & create them as needed"
[ -d $OUT_DIR ] || mkdir -p $OUT_DIR
[ -d $TMP_DIR ] || mkdir -p $TMP_DIR

echo "### Step 0: Index reference" # Run only once! Skip this step if you have already generated BWA index files
bwa index -a bwtsw $REF

echo "### Step 1.A: FASTQ to BAM (1st)"
bwa mem -t $CPU $REF $IN_DIR/DTG-HiC-402_R1_001.fastq.gz | samtools view -@ $CPU -Sb - > $OUT_DIR/DTG-HiC-402_R1_001.bam

echo "### Step 1.B: FASTQ to BAM (2nd)"
bwa mem -t $CPU $REF $IN_DIR/DTG-HiC-402_R2_001.fastq.gz | samtools view -@ $CPU -Sb - > $OUT_DIR/DTG-HiC-402_R2_001.bam

echo "### Step 2.A: Filter 5' end (1st)"
samtools view -h $OUT_DIR/DTG-HiC-402_R1_001.bam | perl $FILTER | samtools view -Sb - > $OUT_DIR/DTG-HiC-402_R1_001_filt.bam

echo "### Step 2.B: Filter 5' end (2nd)"
samtools view -h $OUT_DIR/DTG-HiC-402_R2_001.bam | perl $FILTER | samtools view -Sb - > $OUT_DIR/DTG-HiC-402_R2_001_filt.bam

echo "### Step 3A: Pair reads & mapping quality filter"
perl $COMBINER $OUT_DIR/DTG-HiC-402_R1_001_filt.bam $OUT_DIR/DTG-HiC-402_R2_001_filt.bam samtools $MAPQ_FILTER | samtools view -bS -t $FAIDX - | samtools sort -@ $CPU -o $TMP_DIR/DTG-HiC-402_tmp.bam -

echo "### Step 3.B: Add read group"
java -Xmx4G -Djava.io.tmpdir=temp/ -jar $PICARD AddOrReplaceReadGroups INPUT=$TMP_DIR/DTG-HiC-402_tmp.bam OUTPUT=$OUT_DIR/DTG-HiC-402_paired.bam ID=402 LB=402 SM=HiC PL=ILLUMINA PU=none

echo "### Step 4: Mark duplicates"
java -Xmx30G -XX:-UseGCOverheadLimit -Djava.io.tmpdir=temp/ -jar $PICARD MarkDuplicates INPUT=$OUT_DIR/DTG-HiC-402_paired.bam OUTPUT=$OUT_DIR/$REP_LABEL.bam METRICS_FILE=$OUT_DIR/metrics.$REP_LABEL.txt TMP_DIR=$TMP_DIR ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=TRUE

cd $OUT_DIR

samtools index $OUT_DIR/$REP_LABEL.bam

perl $STATS $OUT_DIR/$REP_LABEL.bam > $OUT_DIR/$REP_LABEL.bam.stats

echo "Finished Mapping Pipeline through Duplicate Removal"

#########################################################################################################################################
###                                       How to Accommodate Biological Replicates                                                    ###
### This pipeline is currently built for processing a single sample with one read1 and read2 fastq file.                              ###
### Biological replicates (eg. multiple libraries made from the same sample) should be merged before proceeding with subsequent steps.###
### The code below is an example of how to merge biological replicates.                                                               ###
#########################################################################################################################################
#
#	INPUTS_BIOLOGICAL_REPS=('bash' 'array' 'of' 'bams' 'from' 'replicates') #BAM files you want combined as biological replicates
#   example bash array - INPUTS_BIOLOGICAL_REPS=('INPUT=A_rep1.bam' 'INPUT=A_rep2.bam' 'INPUT=A_rep3.bam')
#
#	java -Xmx8G -Djava.io.tmpdir=temp/ -jar $PICARD MergeSamFiles $INPUTS_BIOLOGICAL_REPS OUTPUT=$MERGE_DIR/$LABEL.bam USE_THREADING=TRUE ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT
#
#	samtools index $MERGE_DIR/$LABEL.bam

# perl $STATS $MERGE_DIR/$LABEL.bam > $MERGE_DIR/$LABEL.bam.stats

# echo "Finished Mapping Pipeline through merging Biological Replicates"
