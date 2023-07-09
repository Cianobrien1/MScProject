import os
import glob
from rdkit import Chem
from rdkit.Chem import rdMolAlign, AllChem
from collections import defaultdict
from joblib import Parallel, delayed

# Define the directories
conformer_dir = "/home/s2451611/MScProject/pdb_10_percent_alignment_test_sample"
crystal_dir = "/home/s2451611/MScProject/10_percent_crystal_ligand_pdb_db"
target_dir = "/home/s2451611/MScProject/10_aligned_pdbs_db"

# Create the target directory if it doesn't exist
os.makedirs(target_dir, exist_ok=True)

alignment_error_file = "alignment_error.txt"

def align_and_save(pdb_file):
    pdb_id = os.path.basename(pdb_file).split('_')[0]
    mol_conformer = Chem.MolFromPDBFile(pdb_file)

    if mol_conformer is None:
        with open(alignment_error_file, "a") as error_file:
            error_file.write(f"Unable to read conformer: {pdb_file}\n")
        return None

    # Add hydrogens to the conformer
    #mol_conformer = Chem.AddHs(mol_conformer)

    # Read the crystal pose for the ligand
    crystal_file = os.path.join(crystal_dir, pdb_id + "_ligand.pdb")
    mol_crystal = Chem.MolFromPDBFile(crystal_file)

    if mol_crystal is None:
        with open(alignment_error_file, "a") as error_file:
            error_file.write(f"Unable to read crystal: {crystal_file}\n")
        return None

    # Align the conformer to the crystal pose and get the RMSD
    try:
        rmsd = rdMolAlign.AlignMol(mol_conformer, mol_crystal)
    except RuntimeError:
        with open(alignment_error_file, "a") as error_file:
            error_file.write(f"No sub-structure match found between {pdb_file} and {crystal_file}\n")
        return None

    # Write the aligned conformer to a new PDB file
    Chem.MolToPDBFile(mol_conformer, os.path.join(target_dir, os.path.basename(pdb_file)))
    return rmsd

# List of pdb files
pdb_files = glob.glob(conformer_dir + "/*.pdb")

# Run the alignment in parallel using joblib
rmsd_values = Parallel(n_jobs=-1)(delayed(align_and_save)(pdb_file) for pdb_file in pdb_files)

# Filter out the None values
rmsd_values = [x for x in rmsd_values if x is not None]

# Calculate the average RMSD across all conformers
average_rmsd = sum(rmsd_values) / len(rmsd_values)

# Write the average RMSD to a text file
with open("rdkit_rmsd.txt", "w") as f:
    f.write(f"Average RMSD: {average_rmsd}\n")