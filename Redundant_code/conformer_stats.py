##############################################################################################
# This script counts how many conformers have been created for each ligand and saves it to a #
# txt file called "conformer_stats.txt"                                                      #
##############################################################################################

import os
from collections import defaultdict
from rdkit import Chem

# Define the path to your directory
dir_path = "/home/s2451611/MScProject/rdkit_2023_UFF_conformer_dir"

# Create a dictionary to store the number of files for each number of conformers
conformer_dict = defaultdict(int)

# Create a dictionary to store the number of conformers for each PDB code
pdb_dict = {}

# Initialize a variable to store total number of conformers across all ligands
total_conformers = 0

# Loop over each file in the directory
for file_name in os.listdir(dir_path):
    if file_name.endswith('.sdf'):
        pdb_code = file_name[:4]
        file_path = os.path.join(dir_path, file_name)
        
        # Read the .sdf file
        supplier = Chem.SDMolSupplier(file_path)
        mols = [mol for mol in supplier if mol is not None]
        
        # Store the number of conformers for this PDB code
        pdb_dict[pdb_code] = len(mols)
        
        # Update the count of files with this number of conformers
        conformer_dict[len(mols)] += 1

        # Update total conformers
        total_conformers += len(mols)

# Write the results to a file
with open('2023_UFF_conformer_stats.txt', 'w') as f:
    # Write total number of conformers across all ligands at the top
    f.write(f"Total number of conformers across all ligands: {total_conformers}\n\n")

    for num_conformers, num_files in conformer_dict.items():
        f.write(f"{num_files} ligands have {num_conformers} conformers\n")

    f.write("\n")

    for pdb_code, num_conformers in pdb_dict.items():
        f.write(f"{pdb_code} has {num_conformers} conformers\n")
