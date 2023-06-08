from rdkit import Chem
import os

molecule = Chem.rdmolfiles.MolFromPDBFile("/home/s2451611/MScProject/SMILE_pdbs/5ey4_SMILE.pdb", sanitize=False)
Chem.rdmolfiles.MolToMolFile(molecule, "test.mol")