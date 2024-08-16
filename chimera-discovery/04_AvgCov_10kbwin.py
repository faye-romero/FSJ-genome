## Python script that calculates average coverage (1 kb windows) from many samtools mpileup outputs ##
## Usage: run this within a directory with many truncated mpileup output .txt files ##
## 29 July 2024 ##

import pandas as pd
import glob
import os

# Define the file names of the mpileup files you want to consider (I created truncated mpileup outputs using truncate_mpileup_output.py)
file_pattern = 'trunc_H*_aln_sorted_filt_dedup_two_mismatch_cov_Wunloc3_AllPos.txt'
input_directory = '.' # Set wd

# List of all files matching the pattern
file_list = glob.glob(os.path.join(input_directory, file_pattern))

# Set window size
window_size = 1000

for input_file in file_list:
    # Get the sample ID from the filename
    filename = os.path.basename(input_file)
    sample_id = filename.split('_')[1]  # Assuming the sample ID is the second component when split by '_'
    # Read the input truncated mpileup file
    data = pd.read_csv(input_file, delimiter='\t', header=None, names=['Chromosome', 'Position', 'Coverage'])
    output_data = []
    # Group by chromosome to handle multiple chromosomes separately
    grouped = data.groupby('Chromosome')
    for chromosome, group in grouped:
        start_pos = group['Position'].min()
        end_pos = start_pos + window_size
        while start_pos <= group['Position'].max():
            window_data = group[(group['Position'] >= start_pos) & (group['Position'] < end_pos)]
            if not window_data.empty:
                avg_coverage = window_data['Coverage'].mean()
                output_data.append([chromosome, start_pos, end_pos-1, avg_coverage])
            start_pos = end_pos
            end_pos += window_size
    # Convert the output data to a DataFrame
    output_df = pd.DataFrame(output_data, columns=['Chromosome', 'StartPos', 'EndPos', 'AverageCoverage'])
    output_file = os.path.join(input_directory, f'{sample_id}_AvgCov_Wunloc3_AllPos_1kbwin.txt')
    # Write
    output_df.to_csv(output_file, sep='\t', index=False)