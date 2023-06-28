################################################################################## 
# This script copies the ligands from the pdb_crystal poses that have a matching #
# pdb id to the 10_percent_alignment_test_sample                                 #
##################################################################################

import os
import shutil

# Define directories
src_dir = "/home/s2451611/MScProject/openbabel_pdb_crystal_pose_ligands"
dest_dir = "/home/s2451611/MScProject/10_percent_crystal_ligand_pdb"
test_sample_dir = "/home/s2451611/MScProject/10_percent_alignment_test_sample"

# Create destination directory if it doesn't exist
if not os.path.exists(dest_dir):
    os.makedirs(dest_dir)

# Get list of file names in the test_sample_dir
test_sample_files = os.listdir(test_sample_dir)
test_sample_ids = [file.split('_')[0] for file in test_sample_files if file.endswith('.sdf')]

# Iterate over the files in source directory
for file_name in os.listdir(src_dir):
    # Check if the file ends with .sdf
    if file_name.endswith('.pdb'):
        # Get the pdb_id from the file_name
        pdb_id = file_name.split('_')[0]
        # If the pdb_id exists in the test_sample_ids, copy the file
        if pdb_id in test_sample_ids:
            # Define source and destination paths
            src_file_path = os.path.join(src_dir, file_name)
            dest_file_path = os.path.join(dest_dir, file_name)
            
            # Copy the file
            shutil.copy2(src_file_path, dest_file_path)

print("File copying completed.")
