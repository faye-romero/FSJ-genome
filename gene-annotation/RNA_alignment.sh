#!/bin/bash
#SBATCH -p rosalind
#SBATCH --job-name="RNA.alignment"
#SBATCH -N 1
#SBATCH -c 21
#SBATCH --mem=50G
#SBATCH --time=3-00:00:00
#SBATCH --output=RNA.alignment.03152024.log

## RNA read alignment with STAR ##
## Faye Romero. University of Rochester. 18 March 2024 ##

# Load modules
module load rnastar/2.7.3a
module load picard/2.27.4
module load samtools/1.9

# Set wd and variables
cd /scratch/nchen11_lab/FSJ_GeneAnnotation

HOME='/scratch/nchen11_lab/FSJ_GeneAnnotation'
REF='AphCoe_V3_internal_Mar2024.fasta.masked'
RNASEQ='/scratch/nchen11_lab/RNAseq'

## Step 1. Align raw RNA-seq fq files to SOFTMASKED reference genome using STAR ##
#Note: Raw RNA-seq reads have already been processed through TrimGalore.

mkdir -p FSJ_indices_final

#Step 1a. Generate genome indices
echo "generating genome indices. $(date)."
echo "STAR --runMode genomeGenerate \
     --genomeDir ${HOME}/FSJ_indices_final \
     --genomeFastaFiles ${REF} \
     --runThreadN 21"
     
STAR --runMode genomeGenerate \
     --genomeDir ${HOME}/FSJ_indices_final \
     --genomeFastaFiles ${REF} \
     --runThreadN 21

echo "."
echo "."
echo "."

#Step 1b. Align reads

mkdir -p star_alns_March2024

echo "Aligning raw RNA-seq .fq files to the softmasked reference. Using --outSAMstrandField and intronMotif. $(date)."
echo ""
echo "for sample in 243_liver_ACTTGA_L003 243_liver_ACTTGA_L005 251_heart_GATCAG_L003 251_heart_GATCAG_L005 251_kidney_TTAGGC_L003 251_kidney_TTAGGC_L005 253_heart_GCCAAT_L003 253_heart_GCCAAT_L005 253_kidney_ACAGTG_L003 253_kidney_ACAGTG_L005 253_ovary_CGATGT_L003 253_ovary_CGATGT_L005 
do
    STAR --genomeDir ${HOME}/FSJ_indices_final \
        --readFilesIn ${RNASEQ}/NC_FSJ_${sample}_R1_val_1.fq.gz ${RNASEQ}/NC_FSJ_${sample}_R2_val_2.fq.gz \
        --twopassMode Basic --outFileNamePrefix ${HOME}/star_alns_March2024/${sample} \
        --runThreadN 21 --readFilesCommand zcat --outSAMstrandField intronMotif

    java -Xmx25g -jar /software/picard/2.27.4/picard.jar AddOrReplaceReadGroups -VALIDATION_STRINGENCY LENIENT \
        -I ${HOME}/star_alns_March2024/${sample}Aligned.out.sam \
        -O ${HOME}/star_alns_March2024/${sample}.tmp.add.bam \
        -RGID ${sample} -RGLB ${sample} -RGPL RNA -RGPU unit1 -RGSM ${sample}

    samtools sort -n ${HOME}/star_alns_March2024/${sample}.tmp.add.bam > ${HOME}/star_alns_March2024/${sample}.sort.bam
done"
echo "."
echo "."
echo "."

for sample in 243_liver_ACTTGA_L003 243_liver_ACTTGA_L005 251_heart_GATCAG_L003 251_heart_GATCAG_L005 251_kidney_TTAGGC_L003 251_kidney_TTAGGC_L005 253_heart_GCCAAT_L003 253_heart_GCCAAT_L005 253_kidney_ACAGTG_L003 253_kidney_ACAGTG_L005 253_ovary_CGATGT_L003 253_ovary_CGATGT_L005 
do
	echo "aligning ${sample}. $(date)."
    STAR --genomeDir ${HOME}/FSJ_indices_final \
        --readFilesIn ${RNASEQ}/NC_FSJ_${sample}_R1_val_1.fq.gz ${RNASEQ}/NC_FSJ_${sample}_R2_val_2.fq.gz \
        --twopassMode Basic --outFileNamePrefix ${HOME}/star_alns_March2024/${sample} \
        --runThreadN 21 --readFilesCommand zcat --outSAMstrandField intronMotif
	echo ""
	echo "assigning read groups for ${sample}. $(date)."
    java -Xmx25g -jar /software/picard/2.27.4/picard.jar AddOrReplaceReadGroups -VALIDATION_STRINGENCY LENIENT \
        -I ${HOME}/star_alns_March2024/${sample}Aligned.out.sam \
        -O ${HOME}/star_alns_March2024/${sample}.tmp.add.bam \
        -RGID ${sample} -RGLB ${sample} -RGPL RNA -RGPU unit1 -RGSM ${sample}
	echo ""
	echo "samtools sorting for ${sample}. $(date)."
    samtools sort -n ${HOME}/star_alns_March2024/${sample}.tmp.add.bam > ${HOME}/star_alns_March2024/${sample}.sort.bam
    echo "done processing ${sample}."
    echo "."
    echo "."
    echo "."
done