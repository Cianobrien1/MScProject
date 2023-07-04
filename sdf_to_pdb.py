#########################################
# This script converts sdf to pdb files #
#########################################

import os
from rdkit import Chem
from rdkit.Chem import AllChem

# Define the source and target directories
src_dir = "/home/s2451611/MScProject/Raw_Data/PDBBind_ligand_sdf"
tgt_dir = "/home/s2451611/MScProject/Raw_Data/PDBBind_ligand_pdb"

# Create the target directory if it doesn't exist
os.makedirs(tgt_dir, exist_ok=True)

# Iterate over the .sdf files in the source directory
for filename in os.listdir(src_dir):
    if filename.endswith(".sdf"):
        # Determine the base name (without extension) for the file
        basename = os.path.splitext(filename)[0]

        # Read the multi-conformer molecule from the .sdf file
        supplier = Chem.SDMolSupplier(os.path.join(src_dir, filename), sanitize=False)

        # Iterate over the conformers in the molecule
        for i, mol in enumerate(supplier):
            if mol is not None:
                # Generate a filename for the conformer
                pdb_filename = f"{basename}.pdb"

                # Write the conformer to a .pdb file in the target directory
                Chem.MolToPDBFile(mol, os.path.join(tgt_dir, pdb_filename))
