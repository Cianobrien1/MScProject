#########################################
# This script converts sdf to pdb files #
#########################################

import os
from rdkit import Chem
from rdkit.Chem import AllChem

import os
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdchem import KekulizeException

src_dir = "/home/s2451611/MScProject/MOAD_missing_conformers"
tgt_dir = "/home/s2451611/MScProject/MOAD_missing_conformers_pdb"

os.makedirs(tgt_dir, exist_ok=True)

for filename in os.listdir(src_dir):
    if filename.endswith(".sdf"):
        basename = os.path.splitext(filename)[0]
        supplier = Chem.SDMolSupplier(os.path.join(src_dir, filename), sanitize=False)
        for i, mol in enumerate(supplier):
            if mol is not None:
                try:
                    pdb_filename = f"{basename}.pdb"
                    Chem.MolToPDBFile(mol, os.path.join(tgt_dir, pdb_filename))
                except KekulizeException as e:
                    with open("pdb_to_sdf_error.txt", "w") as error_file:
                        error_file.write(f"Could not convert molecule {basename} due to: {e}")

