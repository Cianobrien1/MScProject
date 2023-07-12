##################################################################################################
# This script uses an MCS approach and rdkit to align the 10 percent test samples                #
# to the 10% crystal ligand samples using sdf file inputs for both ligands being aligned         #
##################################################################################################


from rdkit import Chem
from rdkit.Chem import AllChem, rdFMCS
from rdkit.Chem import rdmolops  # import the rdmolops module
import os
import numpy as np
from tqdm import tqdm  # import the tqdm library

# Path to directories
test_sample_dir = "/home/s2451611/MScProject/10_percent_alignment_test_sample"
crystal_dir = "/home/s2451611/MScProject/10_percent_crystal_ligand_sdf_db"
aligned_dir = "/home/s2451611/MScProject/MCS_aligned_pdbs"
rmsd_file_path = "/home/s2451611/MScProject/MCS_rdkit_conformer_rmsd.txt"

# Create directory for aligned pdbs if it doesn't exist
os.makedirs(aligned_dir, exist_ok=True)

rmsds = []

# Open file to write rmsd values
with open(rmsd_file_path, 'w') as rmsd_file:
    for filename in tqdm(os.listdir(test_sample_dir)):  # wrap os.listdir with tqdm
        try:  # add a try block
            pdb_id = filename.split('_')[0]  # Extract pdb_id from filename

            # Load the crystal ligand
            crystal = Chem.SDMolSupplier(os.path.join(crystal_dir, f"{pdb_id}_ligand.sdf"), sanitize=False)[0]
            rdmolops.SanitizeMol(crystal)  # sanitize the crystal ligand

            # Load the multiple conformer molecule
            molsupplier = Chem.SDMolSupplier(os.path.join(test_sample_dir, filename), sanitize=False)
            
            min_rmsd = None
            min_conf_id = None
            for i, mol in enumerate(molsupplier):
                rdmolops.SanitizeMol(mol)  # sanitize the conformer
                # Ensure a consistent atom ordering
                mol = Chem.RenumberAtoms(mol, Chem.CanonicalRankAtoms(mol))
                crystal = Chem.RenumberAtoms(crystal, Chem.CanonicalRankAtoms(crystal))

                # Find the MCS between the crystal and the conformer
                mcs = rdFMCS.FindMCS([crystal, mol])
                patt = Chem.MolFromSmarts(mcs.smartsString)
                atom_map = list(zip(crystal.GetSubstructMatch(patt), mol.GetSubstructMatch(patt)))

                try:
                    rmsd = AllChem.AlignMol(mol, crystal, atomMap=atom_map)
                except RuntimeError:
                    print(f"Error aligning conformer {pdb_id}_{i}")
                    continue
                rmsds.append(rmsd)
                rmsd_file.write(f'{pdb_id}_{i}: {rmsd}\n')  # Write the rmsd to the file
                if min_rmsd is None or rmsd < min_rmsd:
                    min_rmsd = rmsd
                    min_conf_id = i
                    
            # Save the conformer with the lowest rmsd
            if min_conf_id is not None:
                mol = molsupplier[min_conf_id]
                Chem.MolToPDBFile(mol, os.path.join(aligned_dir, f"{pdb_id}_aligned.pdb"))
        except Exception as e:  # add an except block
            print(f"Error processing file {filename}: {e}")
            continue

# Write the average rmsd at the top of the file
if rmsds:
    average_rmsd = np.mean(rmsds)
    with open(rmsd_file_path, 'r+') as rmsd_file:
        content = rmsd_file.read()
        rmsd_file.seek(0, 0)
        rmsd_file.write(f"Average RMSD: {average_rmsd}\n" + content)
