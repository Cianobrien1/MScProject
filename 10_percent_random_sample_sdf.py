import os
import random
import re

source_dir = '/home/s2451611/MScProject/conformer_dir' # replace with your source directory
output_file = '10_percent_sample_for_alignment.txt' # replace with your output file path

# Get a list of all SDF files in the source directory
sdf_files = [f for f in os.listdir(source_dir) if f.endswith('.sdf')]

# Calculate 10% of the total number of files
num_test_files = len(sdf_files) // 10

# Randomly select a subset of files
test_files = random.sample(sdf_files, num_test_files)

# Extract PDB IDs and write to the output file
with open(output_file, 'w') as f:
    for file in test_files:
        pdb_id = re.match(r"(\w+)_SMILE.sdf", file).group(1)
        f.write(pdb_id + '\n')
