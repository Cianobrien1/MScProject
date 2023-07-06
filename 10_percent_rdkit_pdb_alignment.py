##################################################################################################
# This script uses rdkit to align the 10 percent test samples to the 10% crystal ligand samples  #
# using pdb file inputs for both ligands being aligned                                           #
# Hydrogens are added only to the conformer                                                      #
##################################################################################################

import os
import glob
from rdkit import Chem
from rdkit.Chem import rdMolAlign, AllChem
from collections import defaultdict
from joblib import Parallel, delayed

# Define the directories
conformer_dir = "/home/s2451611/MScProject/MOAD_missing_conformers_pdb"
crystal_dir = "/home/s2451611/MScProject/MOAD_missing_pdbs"
target_dir = "/home/s2451611/MScProject/missing_aligned_pdbs"

# Create the target directory if it doesn't exist
os.makedirs(target_dir, exist_ok=True)

def align_and_save(pdb_file):
    pdb_id = os.path.basename(pdb_file).split('_')[0]
    mol_conformer = Chem.MolFromPDBFile(pdb_file)

    # Add hydrogens to the conformer
    mol_conformer = Chem.AddHs(mol_conformer)

    # Read the crystal pose for the ligand
    crystal_file = os.path.join(crystal_dir, pdb_id + "_ligand.pdb")
    mol_crystal = Chem.MolFromPDBFile(crystal_file)

    # Align the conformer to the crystal pose and get the RMSD
    rmsd = rdMolAlign.AlignMol(mol_conformer, mol_crystal)

    # Write the aligned conformer to a new PDB file
    Chem.MolToPDBFile(mol_conformer, os.path.join(target_dir, os.path.basename(pdb_file)))
    return rmsd

# List of pdb files
pdb_files = glob.glob(conformer_dir + "/*.pdb")

# Run the alignment in parallel using joblib
rmsd_values = Parallel(n_jobs=-1)(delayed(align_and_save)(pdb_file) for pdb_file in pdb_files)

# Calculate the average RMSD across all conformers
average_rmsd = sum(rmsd_values) / len(rmsd_values)

# Write the average RMSD to a text file
with open("rdkit_rmsd.txt", "w") as f:
    f.write(f"Average RMSD: {average_rmsd}\n")
