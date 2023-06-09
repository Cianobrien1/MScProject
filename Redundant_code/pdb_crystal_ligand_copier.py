#############################################################################################
# This script was used to copy all the ligand pdb files from the "pdb_all_receptor_ligands" #
# directory to thier own directory                                                          #
#############################################################################################

import os
import glob
import shutil

output_dir = "/home/s2451611/MScProject/pdb_all_receptor_ligands"
new_dir = "/home/s2451611/MScProject/crystal_pose_ligands_pdb"

# Get a list of all PDB files with 'ligand' in the name in the output directory and its subdirectories
pdb_files = glob.glob(os.path.join(output_dir, "**/*ligand*.pdb"), recursive=True)

print(f"Found {len(pdb_files)} ligand PDB files.")

# Create the new directory if it doesn't exist
os.makedirs(new_dir, exist_ok=True)

# Copythe files
for pdb_file in pdb_files:
    # Get the base name of the file
    base_name = os.path.basename(pdb_file)
    # Construct the new file path
    new_file_path = os.path.join(new_dir, base_name)
    # Copy the file
    shutil.copy(pdb_file, new_file_path)

print(f"Copied all ligand PDB files to {new_dir}.")
