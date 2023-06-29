##################################################################################################
# This script uses rdkit to align the 10 percent test samples to the 10% crystal ligand samples  #
# using pdb file inputs for both ligands being aligned                                           #
##################################################################################################

import os
import glob
from rdkit import Chem
from rdkit.Chem import rdMolAlign
from collections import defaultdict

# Define the directories
conformer_dir = "/home/s2451611/MScProject/pdb_10_percent_alignment_test_sample"
crystal_dir = "/home/s2451611/MScProject/10_percent_crystal_ligand_pdb"
target_dir = "/home/s2451611/MScProject/10_percent_aligned_pdb"

# Create the target directory if it doesn't exist
os.makedirs(target_dir, exist_ok=True)

# Initialize a dictionary to store the RMSD values
rmsd_dict = defaultdict(list)

# Iterate over the conformer PDB files
for pdb_file in glob.glob(conformer_dir + "/*.pdb"):
    pdb_id = os.path.basename(pdb_file).split('_')[0]
    mol_conformer = Chem.MolFromPDBFile(pdb_file)

    # Read the crystal pose for the ligand
    crystal_file = os.path.join(crystal_dir, pdb_id + "_ligand.pdb")
    mol_crystal = Chem.MolFromPDBFile(crystal_file)

    # Align the conformer to the crystal pose and get the RMSD
    rmsd = rdMolAlign.AlignMol(mol_conformer, mol_crystal)
    rmsd_dict[pdb_id].append((rmsd, pdb_file))

# Write the aligned conformer with the lowest RMSD for each ligand to a new PDB file
for pdb_id, rmsd_list in rmsd_dict.items():
    min_rmsd, min_rmsd_file = min(rmsd_list)
    mol = Chem.MolFromPDBFile(min_rmsd_file)
    Chem.MolToPDBFile(mol, os.path.join(target_dir, os.path.basename(min_rmsd_file)))

# Calculate the average RMSD across all conformers for all ligands
total_rmsd = sum(rmsd for rmsd_list in rmsd_dict.values() for rmsd, _ in rmsd_list)
total_conformers = sum(len(rmsd_list) for rmsd_list in rmsd_dict.values())
average_rmsd = total_rmsd / total_conformers

# Write the average RMSD to a text file
with open("rdkit_rmsd.txt", "w") as f:
    f.write(f"Average RMSD: {average_rmsd}\n")
