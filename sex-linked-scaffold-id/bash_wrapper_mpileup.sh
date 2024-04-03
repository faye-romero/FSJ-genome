#!/bin/bash

# Path to the directory containing the subdirectories
BASE_DIR="/scratch/nchen11_lab/new_FSJgenome/sex_linked_scaffs"

# Path to the list of subdirectories
SAMPLE_LIST="list_all_samples.txt"

# Path to the script to be executed
SCRIPT_PATH="/scratch/nchen11_lab/new_FSJgenome/sex_linked_scaffs/mpileup.sh"

# Read the list of subdirectories from the file
while IFS= read -r SAMPLE; do
  # Form the subdirectory name by appending "_out"
  SUBDIR="${SAMPLE}_out"
  SUBDIR_PATH="$BASE_DIR/$SUBDIR"

  # Check if the subdirectory exists
  if [ -d "$SUBDIR_PATH" ]; then
    # Submit a job for each subdirectory using sbatch
    sbatch -J "${SAMPLE}_job" -o "${SUBDIR_PATH}/${SAMPLE}_mpileup.out.txt" -e "${SUBDIR_PATH}/${SAMPLE}_mpileup.err.txt" "$SCRIPT_PATH" "$SUBDIR_PATH" "$SAMPLE"
  else
    echo "Error: Subdirectory $SUBDIR_PATH does not exist."
  fi
done < "$SAMPLE_LIST"

echo "All jobs submitted."
