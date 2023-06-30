#####################################################################################
# This script was used to check rdkit's ability to parse the pdb files being passed #
# to it on an individual basis, and convert them to mol files as confgen originally #
# took mol files as input, not SMILEs as in the custom modified version             #
#####################################################################################

from rdkit import Chem
import os

molecule = Chem.rdmolfiles.MolFromPDBFile("/home/s2451611/MScProject/SMILE_pdbs/5ey4_SMILE.pdb", sanitize=False)
Chem.rdmolfiles.MolToMolFile(molecule, "test.mol")