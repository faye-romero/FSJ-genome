## Python script to truncate the output of samtools mpileup to just the scaffold, position, and coverage ##
## Usage: run this within a directory with many mpileup output .txt files ##
## 29 July 2024 ##

import os
import pandas as pd

# Set wd
directory = '.'

# List all files
files = os.listdir(directory)

for filename in files:
    if filename.endswith('.txt'):
        filepath = os.path.join(directory, filename)
        df = pd.read_csv(filepath, sep='\t', header=None)
        # Extract the first, second, and fourth columns
        truncated_df = df[[0, 1, 3]]
        # Define the new filename with "trunc_" prefix
        new_filename = 'trunc_' + filename
        new_filepath = os.path.join(directory, new_filename)
        # Write
        truncated_df.to_csv(new_filepath, sep='\t', index=False, header=False)