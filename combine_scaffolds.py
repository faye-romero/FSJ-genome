## combine_scaffolds.py: Python script to join multiple scaffolds from a FASTA file in one single large scaffold, each chunk separated by 100 Ns ##
## Faye Romero ##
## 11 Jan 2024 ##

import sys

if len(sys.argv) != 4:
    print("Usage: python combine_scaffolds.py input_file output_file large_scaffold_header")
    print("Example usage: python combine_scaffolds.py scaffolds.fasta combined_scaffolds.fasta '>large_scaffold_1'")
    print("This script takes in a single-line FASTA file with multiple scaffolds and combines the scaffold sequences together under one new scaffold (specified by the argument 'large_scaffold_header'). Each sequence will be separated by 100 Ns.")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]
large_scaffold_header = sys.argv[3]

# Read the input file
with open(input_file, "r") as infile:
    lines = infile.readlines()

# Initialize variables
large_scaffold_sequence = ""

# Process the input data
for line in lines:
    line = line.strip()
    if line.startswith(">"):
        # Skip scaffold headers
        if large_scaffold_sequence:
            large_scaffold_sequence += "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
        continue
    large_scaffold_sequence += line

# Write the output to a new file
with open(output_file, "w") as outfile:
    outfile.write(f"{large_scaffold_header}\n{large_scaffold_sequence}\n")

print(f"Large scaffold created and saved to {output_file}.")
