########################################################################################################
# This script moves all pdb files for ligands in the active dataset from the Iridium dataset           #
# directory to the new directory /home/s2451611/MScProject/Raw_Data/Active_ligand_crystal_poses_pdb_db #
########################################################################################################

import os
import shutil

# Define paths
input_filepath = '/home/s2451611/MScProject/Data_info/Data_Source/Iridium.txt'
data_directory = '/home/s2451611/MScProject/Raw_Data/Iridium_Successes'
output_directory = '/home/s2451611/MScProject/Raw_Data/Iridum_only_ligand_pdbs'
missing_filepath = '/home/s2451611/MScProject/Iridium_missing_pdb_ids.txt'

# Check and create the output directory if it does not exist
if not os.path.exists(output_directory):
    os.makedirs(output_directory)

# Read pdb ids from the file
with open(input_filepath, 'r') as file:
    pdb_ids = [line.strip() for line in file]

# Open the missing pdb ids file
with open(missing_filepath, 'w') as missing_file:
    # Iterate over all pdb ids
    for pdb_id in pdb_ids:
        # Construct file path
        pdb_file_path = os.path.join(data_directory, pdb_id, f'{pdb_id}_ligand.pdb')

        # Check if file exists
        if os.path.isfile(pdb_file_path):
            # Construct output file path
            output_file_path = os.path.join(output_directory, f'{pdb_id}_ligand.pdb')

            # Copy file
            shutil.copy2(pdb_file_path, output_file_path)
        else:
            # Write missing pdb id to the file
            missing_file.write(pdb_id + '\n')
