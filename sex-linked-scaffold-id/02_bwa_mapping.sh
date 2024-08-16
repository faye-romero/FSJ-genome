#!/bin/bash
#SBATCH -p rosalind
#SBATCH --job-name="bwa.mapping"
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=30G
#SBATCH --nodelist=bhd0057
#SBATCH --time=7-00:00:00
#SBATCH --output=bwamem.out.log
#SBATCH --error=bwamem.err.log

## PART 2 SLURM script for finding sex-linked scaffolds given a female reference genome and one set of high-coverage paired-end Illumina data.##
## Faye Romero. University of Rochester. #
## 13 Feb 2024 ##

# Subdirectory path passed as an argument
SUBDIR_PATH="$1"

# Extract the sample name from the subdirectory path
NAME=$(basename "$SUBDIR_PATH" | cut -d '_' -f 1)

# Set absolute working directory
wd=/scratch/nchen11_lab/Faye/new_FSJgenome/sex_linked_scaffs/${NAME}_out
cd $wd

# Let's begin!
echo "Processing sample ${NAME} within the working dir ${wd}. Let's gooooo!"
echo "$(date)."
echo "................................."

################################
####### Step 2. BWA mem ########
################################

#SBATCH --job-name="bwamem"
#SBATCH --output=bwamem.out.log
#SBATCH --error=bwamem.err.log

# Load modules
module load bwa/0.7.4
module load samtools/1.9

# Make sure you're in the proper working dir
cd $wd

# Set variables
REF='/scratch/nchen11_lab/Faye/new_FSJgenome/FSJgenome_July2024_FINAL.fasta'
trim_out_dir="${wd}/trimgalore_out"

# Create output directory, if it does not exist
mkdir -p "${wd}/bwa_out"
bwa_out_dir="${wd}/bwa_out"

# Identify read pairs in trimgalore_out
		pairs=(./trimgalore_out/*_val_1.fq)
        for trimfq1 in "${pairs[@]}"; do
            trimfq2="${trimfq1/_1_val_1/_2_val_2}"
			# Set output file names
			prefix=$(echo "$trimfq1" | cut -d'/' -f3 | sed 's/_1_val_1.fq//')
			echo "Output file prefix will be:"
			echo "$prefix"
			echo ".............................................." 
			
			# Extract read group info from one of the fastq files
			header=$(zcat "$trimfq1" | head -n 1) #prints fq header
			id=$(echo "$header" | cut -f 3-4 -d":" | sed 's/@//' | sed 's/:/./g') #prints flow cell ID and flow cell lane
			sm=$(echo "$trimfq1" | cut -d'/' -f3 | cut -d'_' -f1) #extracts sample ID from fq
			lb="${sm}_lib" #prints sample ID + _lib suffix
			rg=$(echo "@RG\tID:${id}\tSM:${sm}\tLB:${lb}\tPL:ILLUMINA")
			echo "Read group input for $prefix:"
			echo "$rg"
			
			# Run bwa mem for trimfq1 and trimfq2
			echo ".....Running bwa mem on $trimfq1 and $trimfq2....."
			echo "$(date).....bwa mem -M -t 4 -R $(echo "@RG\tID:${id}\tSM:${sm}\tLB:${lb}\tPL:ILLUMINA") "$REF" "$trimfq1" "$trimfq2" > ${prefix}.sam 2> ${prefix}.bwa.out....."
		        bwa mem -M -t 4 -R $(echo "@RG\tID:${id}\tSM:${sm}\tLB:${lb}\tPL:ILLUMINA") "$REF" "$trimfq1" "$trimfq2" > ${prefix}.sam 2> ${prefix}.bwa.out

			# Sort bam file
			echo ".....Sort bam file....."
			echo "$(date).....samtools sort -o ${prefix}.sorted.bam ${prefix}.sam ....."
			samtools sort -o ${prefix}.sorted.bam ${prefix}.sam
			
			# Delete sam file
			rm -f ${prefix}.sam
			
        done

echo "$(date).....bwa mem has ended....."
echo ".............................................."
echo "$(date).....merge and index bam files....."

echo "samtools merge ${NAME}.bam *.bam"
samtools merge ${NAME}.bam *.sorted.bam

echo "samtools index ${NAME}.bam"
samtools index ${NAME}.bam

# move final output to output directory
mv ${NAME}.bam "$bwa_out_dir"
mv ${NAME}.bam.bai "$bwa_out_dir"

# # Check if the bam alignment exists and is not empty
 if [ ! -s "$bwa_out_dir"/"$NAME".bam ] || [ ! -s "$bwa_out_dir"/"$NAME".bam.bai ]; then
     echo "Error: The aligned, sorted bam file and/or the bam index file is missing or empty in $bwa_out_dir!"
     exit 1  # Terminate the job with an error status
 fi

# Delete intermediate files
echo ".............................................."
echo "$(date).....deleting intermediate bam files....."
echo "rm -rf ./*.sorted.bam"
rm -rf ./*.sorted.bam

echo "$(date).....transforming sam to bam, sorting, indexing has ended....."

