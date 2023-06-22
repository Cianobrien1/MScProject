import os
import shutil
import random

source_dir = '/home/s2451611/MScProject/conformer_dir' # replace with your source directory
destination_dir = '/home/s2451611/MScProject/10_percent_alignment_test_sample' # replace with your destination directory

# Get a list of all SDF files in the source directory
sdf_files = [f for f in os.listdir(source_dir) if f.endswith('.sdf')]

# Calculate 10% of the total number of files
num_test_files = len(sdf_files) // 10

# Randomly select a subset of files
test_files = random.sample(sdf_files, num_test_files)

# Copy the selected files to the destination directory
for file in test_files:
    shutil.copy(os.path.join(source_dir, file), destination_dir)
