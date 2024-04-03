#!/bin/bash
#SBATCH -p standard
#SBATCH --job-name="bwa.mapping"
#SBATCH -N 1
#SBATCH -c 12
#SBATCH --mem=40G
#SBATCH --time=5-00:00:00
#SBATCH --output=bwamem.out.log
#SBATCH --error=bwamem.err.log

## PART 2 SLURM script for finding sex-linked scaffolds given a female reference genome and a batch of high-coverage paired-end Illumina data from multiple individuals. ##
## Faye Romero. University of Rochester. #
## 13 Feb 2024 ##

# To run this script, run its wrapper: sh bash_wrapper_bwa_mapping_v2.sh

# Subdirectory path passed as an argument
SUBDIR_PATH="$1"

# Extract the sample name from the subdirectory path
NAME=$(basename "$SUBDIR_PATH" | cut -d '_' -f 1)

# Set absolute working directory
wd=/scratch/nchen11_lab/new_FSJgenome/sex_linked_scaffs/${NAME}_out
cd $wd

# Let's begin!
echo "Processing ${NAME} within the working dir ${wd}. Let's gooooo!"
echo "$(date)."
echo "."
echo "."
echo "."

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
REF='/scratch/nchen11_lab/new_FSJgenome/allmaps_out_v2/FSJv3_internal_Feb2024.fasta'
# export PATH=/scratch/nchen11_lab/new_FSJgenome:$PATH # make sure index files are accessible
trim_out_dir="${wd}/trimgalore_out"

# Create output directory, if it does not exist
mkdir -p "${wd}/bwa_out"
bwa_out_dir="${wd}/bwa_out"

# Make a symlink to the reference index files
cd $bwa_out_dir
ln -s /scratch/nchen11_lab/new_FSJgenome/allmaps_out_v2/FSJv3_internal_Feb2024.fasta.sa
ln -s /scratch/nchen11_lab/new_FSJgenome/allmaps_out_v2/FSJv3_internal_Feb2024.fasta.pac
ln -s /scratch/nchen11_lab/new_FSJgenome/allmaps_out_v2/FSJv3_internal_Feb2024.fasta.ann
ln -s /scratch/nchen11_lab/new_FSJgenome/allmaps_out_v2/FSJv3_internal_Feb2024.fasta.amb
ln -s /scratch/nchen11_lab/new_FSJgenome/allmaps_out_v2/FSJv3_internal_Feb2024.fasta.bwt
ln -s /scratch/nchen11_lab/new_FSJgenome/allmaps_out_v2/FSJv3_internal_Feb2024.fasta
cd $wd

# Identify read pairs in trimgalore_out
        pairs=(trimgalore_out/*_val_1.fq)
        echo ".....Identified pairs of fq files for sample $NAME = $pairs....."
        for pair1 in "${pairs[@]}"; do
            pair2="${pair1/_1_val_1/_2_val_2}"
            # Check if the second pair exists
            if [ -f "$pair2" ]; then
                # Run bwa mem with both pairs
                output="${pair1/_1_val_1.fq/.sam}"
                bamout="${pair1/_1_val_1.fq/.bam}"
                echo ".....Running bwa mem on $pair1 and $pair2....."
                echo "$(date).....bwa mem -t 12 $REF $pair1 $pair2 > $output....."
                bwa mem -t 12 "$REF" "$pair1" "$pair2" > "$output"
                samtools view -@ 12 -Sb -o "$bamout" "$output"
                mv "$output" "$bamout" "$bwa_out_dir"
            else
                # Run bwa mem with single pair
                output="${pair1/_1_val_1.fq/.sam}"
                bamout="${pair1/_1_val_1.fq/.bam}"
                echo ".....Running bwa on $pair1....."
                echo "$(date).....bwa mem -t 12 $REF $pair1 > $output....."
                bwa mem -t 12 "$REF" "$pair1" > "$output"
                samtools view -@ 12 -Sb -o "$bamout" "$output"
                mv "$output" "$bamout" "$bwa_out_dir"
            fi
        done

echo "$(date).....bwa mem has ended....."
echo "$(date).....transforming sam to bam and sorting....."

cd $bwa_out_dir
samtools merge ${NAME}_aln.bam *.bam
# samtools view -@ 12 -Sb -o ${bwa_out_dir}/${NAME}_aln.bam ${bwa_out_dir}/${NAME}_aln.sam
samtools sort -O bam -o ${NAME}_aln_sorted.bam ${NAME}_aln.bam
cd $wd

# Check if the bam alignment exists and is not empty
if [ ! -s "$bwa_out_dir"/"$NAME"_aln_sorted.bam ]; then
    echo "Error: The sorted bam alignment file is missing or empty in $bwa_out_dir!"
    exit 1  # Terminate the job with an error status
fi

# Delete intermediate files
echo "$(date).....deleting intermediate sam and bam files....."
rm -rf ${bwa_out_dir}/${NAME}_aln.sam
rm -rf ${bwa_out_dir}/${NAME}_aln.bam

# Index sorted bam file
echo "$(date).....indexing the sorted bam file....."
cd $bwa_out_dir
samtools index ${NAME}_aln_sorted.bam
cd $wd

echo "$(date).....transforming sam to bam, sorting, indexing has ended....."

# Check if file exist in the "bwa_out" subdirectory and is not empty
if [ ! -s "$bwa_out_dir"/"$NAME"_aln_sorted.bam ] || [ ! -s "$bwa_out_dir"/"$NAME"_aln_sorted.bam.bai ]; then
    echo "Error: The aligned, sorted bam file and/or the BAI index file is missing or empty in $bwa_out_dir!"
fi
