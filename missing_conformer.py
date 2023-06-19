import os

# Define the directories
active_smiles_dir = "/home/s2451611/MScProject/Acitve_Smiles_pdbs"
conformer_dir = "/home/s2451611/MScProject/conformer_dir"

# Get the list of files in each directory
active_smiles_files = os.listdir(active_smiles_dir)
conformer_files = os.listdir(conformer_dir)

# Extract the 4-digit identifiers from the file names in each directory
active_smiles_ids = {f[:4] for f in active_smiles_files}
conformer_ids = {f[:4] for f in conformer_files}

# Find the ids that are in active_smiles_ids but not in conformer_ids
missing_ids = active_smiles_ids - conformer_ids

# Print the missing ids
for id in missing_ids:
    print(id)
