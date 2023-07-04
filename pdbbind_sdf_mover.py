########################################################################################################
# This script moves all sdf files for ligands in the active dataset from the PDBBind dataset directory #
# to the new directory /home/s2451611/MScProject/PDBBind_ligand_sdf/                                   #
########################################################################################################

import os
import shutil

# Define paths
input_filepath = '/home/s2451611/MScProject/Data_info/Data_Source/PDBBind.txt'
data_directory = '/home/s2451611/MScProject/Raw_Data/PDBbind_v2020_refined/refined-set/'
output_directory = '/home/s2451611/MScProject/Raw_Data/PDBBind_ligand_sdf/'
missing_filepath = '/home/s2451611/MScProject/Data_info/Data_Source/missing_pdb_ids.txt'

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
        sdf_file_path = os.path.join(data_directory, pdb_id, f'{pdb_id}_ligand.sdf')

        # Check if file exists
        if os.path.isfile(sdf_file_path):
            # Construct output file path
            output_file_path = os.path.join(output_directory, f'{pdb_id}_ligand.sdf')

            # Copy file
            shutil.copy2(sdf_file_path, output_file_path)
        else:
            # Write missing pdb id to the file
            missing_file.write(pdb_id + '\n')
