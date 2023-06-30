#####################################################################
# This scrit converts all pdbqt files in the "all_receptor_ligands" #
# directory to pdb files using openbabel                            #
#####################################################################

import glob
import os
from openbabel import pybel
from joblib import Parallel, delayed

def convert_pdbqt_to_pdb(file_path, base_dir, output_dir):
    mol = next(pybel.readfile("pdbqt", file_path))
    
    # Get the path to the file, relative to base_dir
    relative_path = os.path.relpath(file_path, base_dir)
    
    # Replace the '.pdbqt' extension with '.pdb'
    output_file_path = os.path.join(output_dir, os.path.splitext(relative_path)[0] + ".pdb")
    
    # Create the output directory if it doesn't exist
    os.makedirs(os.path.dirname(output_file_path), exist_ok=True)
    
    # Write the molecule to a PDB file
    mol.write("pdb", output_file_path, overwrite=True)

base_dir = "/home/s2451611/MScProject/all_receptors_ligands/all_receptors_ligands"
output_dir = "/home/s2451611/MScProject/openbabel_all_receptors_ligands"

# Find all PDBQT files in base_dir and its subdirectories
pdbqt_files = glob.glob(os.path.join(base_dir, "**/*.pdbqt"), recursive=True)

# Run the conversion in parallel using joblib
Parallel(n_jobs=-1)(delayed(convert_pdbqt_to_pdb)(file_path, base_dir, output_dir) for file_path in pdbqt_files)
