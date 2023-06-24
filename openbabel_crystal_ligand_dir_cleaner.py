import os
import glob

conformer_dir = "/home/s2451611/MScProject/conformer_dir"
crystal_dir = "/home/s2451611/MScProject/openbabel_pdb_crystal_pose_ligands"

# Get a list of all SDF files in the conformer_dir
conformer_files = glob.glob(os.path.join(conformer_dir, "*.sdf"))

# Extract PDB IDs from the conformer file names
conformer_pdb_ids = {os.path.basename(f).split('_')[0] for f in conformer_files}

# Get a list of all PDB files in the crystal_dir
crystal_files = glob.glob(os.path.join(crystal_dir, "*.pdb"))

for crystal_file in crystal_files:
    # Extract the PDB ID from the crystal file name
    crystal_pdb_id = os.path.basename(crystal_file).split('_')[0]
    
    # If the crystal PDB ID is not in the set of conformer PDB IDs, remove the file
    if crystal_pdb_id not in conformer_pdb_ids:
        print(f"Removing {crystal_file}...")
        os.remove(crystal_file)

print("Finished removing files.")
