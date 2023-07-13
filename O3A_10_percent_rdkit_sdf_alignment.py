############################################################################################################
# This script uses O3A instead of AllChem.MolAlign approach and rdkit to align the 10 percent test samples #
# to the 10% crystal ligand samples using sdf file inputs for both ligands being aligned                   #
############################################################################################################

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolAlign
from joblib import Parallel, delayed
from tqdm import tqdm
import os
import numpy as np

# Path to directories
test_sample_dir = "/home/s2451611/MScProject/2023_MMFF_10_percent_sdf"
crystal_dir = "/home/s2451611/MScProject/openbabel_10_percent_sdf_crystal_pose"
aligned_dir = "/home/s2451611/MScProject/MMFF_2023_O3A_aligned_pdbs"
rmsd_file_path = "/home/s2451611/MScProject/MMFF_2023_O3A_rdkit_conformer_rmsd.txt"

# Create directory for aligned pdbs if it doesn't exist
os.makedirs(aligned_dir, exist_ok=True)

def process_file(filename):
    rmsds = []

    pdb_id = filename.split('_')[0]  # Extract pdb_id from filename

    # Load the crystal ligand
    crystal = Chem.SDMolSupplier(os.path.join(crystal_dir, f"{pdb_id}_ligand.sdf"))[0]

    # Load the multiple conformer molecule
    molsupplier = Chem.SDMolSupplier(os.path.join(test_sample_dir, filename))
    
    min_rmsd = None
    min_conf_id = None
    for i, mol in enumerate(molsupplier):
        # Ensure a consistent atom ordering
        mol = Chem.RenumberAtoms(mol, Chem.CanonicalRankAtoms(mol))
        crystal = Chem.RenumberAtoms(crystal, Chem.CanonicalRankAtoms(crystal))
        try:
            o3a = rdMolAlign.GetO3A(mol, crystal)  # get the optimal alignment
            rmsd = o3a.Align()  # align the molecules and get the RMSD
        except RuntimeError:
            print(f"Error aligning conformer {pdb_id}_{i}")
            continue
        rmsds.append(rmsd)
        if min_rmsd is None or rmsd < min_rmsd:
            min_rmsd = rmsd
            min_conf_id = i
            
    # Save the conformer with the lowest rmsd
    if min_conf_id is not None:
        mol = molsupplier[min_conf_id]
        Chem.MolToPDBFile(mol, os.path.join(aligned_dir, f"{pdb_id}_aligned.pdb"))
    
    return rmsds

filenames = os.listdir(test_sample_dir)
results = Parallel(n_jobs=-1)(delayed(process_file)(filename) for filename in tqdm(filenames, desc="Processing files"))

all_rmsds = [rmsd for sublist in results for rmsd in sublist]

# Write rmsd values and average rmsd
if all_rmsds:
    with open(rmsd_file_path, 'w') as rmsd_file:
        rmsd_file.write(f"Average RMSD: {np.mean(all_rmsds)}\n")
        for filename, rmsds in zip(filenames, results):
            pdb_id = filename.split('_')[0]
            for i, rmsd in enumerate(rmsds):
                rmsd_file.write(f'{pdb_id}_{i}: {rmsd}\n')
