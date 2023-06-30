#################################################################################
# This script uses openbabel to convert SMILES to pdb files, it was not used as #
# use of openbable here would have mandated re-doing all data with openbabel    #
# instead of rdkit                                                              #
#################################################################################

from openbabel import pybel
import os
from openbabel import openbabel as ob

smile_file = "/home/s2451611/MScProject/Broken_Decoy_Smiles.txt"
i = 0
os.mkdir("openbabel_dir")
with open(smile_file) as file:
    for eachline in file.readlines():
        splitline = eachline.split()
        smile = splitline[6]
        i += 1
        # Create a molecule from the SMILES string
        mol = pybel.readstring("smi", smile)

        # Convert to 3D. This is equivalent to EmbedMolecule in RDKit.
        mol.make3D()

        # If you want to output the molecule in a particular format (e.g., PDB), you can do so:
        pdb_output = mol.write("pdb")
        with open(f"openbabel_dir/{i}.pdb", "w") as f:
            f.write(pdb_output)

    
    
