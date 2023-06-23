from rdkit import Chem
from joblib import Parallel, delayed
import os
import glob

def convert_pdb_to_sdf(pdb_file, sdf_dir):
    # Create an RDKit molecule object from the PDB file
    mol = Chem.rdmolfiles.MolFromPDBFile(pdb_file)
    
    # Construct the output file path
    base_name = os.path.basename(pdb_file).replace('.pdb', '.sdf')
    sdf_file = os.path.join(sdf_dir, base_name)
    
    # Write the molecule to an SDF file
    Chem.rdmolfiles.MolToMolFile(mol, sdf_file)

# Directory containing the input PDB files
pdb_dir = "/home/s2451611/MScProject/pdb_crystal_pose_ligands"

# Directory where the output SDF files will be saved
sdf_dir = "/home/s2451611/MScProject/sdf_crystal_pose_ligands"

# Create the output directory if it doesn't exist
os.makedirs(sdf_dir, exist_ok=True)

# Get a list of all PDB files in the pdb_dir
pdb_files = glob.glob(os.path.join(pdb_dir, "*.pdb"))

# Convert the PDB files to SDF files in parallel
Parallel(n_jobs=-1)(delayed(convert_pdb_to_sdf)(pdb_file, sdf_dir) for pdb_file in pdb_files)
