## extract_unmasked_seq_pos.py: Python script to extract the positions of unmasked (uppercase) regions from a softmasked FASTA file ##
## Faye Romero ##
## 11 Jan 2024 ##

import sys

def process_sequence(scaffold_name, sequence):
    runs = []
    start = None
    for i, char in enumerate(sequence):
        if char.isupper():
            if start is None:
                start = i
        elif start is not None:
            runs.append((start, i - 1))
            start = None
    if start is not None:
        runs.append((start, len(sequence) - 1))

    results = [(scaffold_name, run[0], run[1]) for run in runs]
    return results

def main(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        current_scaffold_name = None
        for line in infile:
            line = line.strip()
            if line.startswith('>'):
                current_scaffold_name = line[1:]
            else:
                results = process_sequence(current_scaffold_name, line)
                for result in results:
                    outfile.write(f"{result[0]} {result[1]} {result[2]}\n")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python process_sequences.py input_file output_file")
        print ("This script processes a softmasked FASTA file (such as the .masked output from RepeatMasker) and identifies runs of consecutive uppercase letters (i.e. unmasked sequence) across the genome. It then outputs the start and end positions of these runs along with the scaffold name. The output is a table (BED format) with the scaffold header, the start position, and the end position of the uppercase runs (0-based). This file is useful as input for samtools mpileup, if you want to ignore softmasked/repetitive regions.")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    main(input_file, output_file)

