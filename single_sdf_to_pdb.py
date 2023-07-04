##################################################################################
# This script converts a single sdf file, including all molecules it may contain #
# # and converts each molecule to a pdb file using rdkit                         #
###################################################################################


from rdkit import Chem
from rdkit.Chem import AllChem

# Read the SDF file
supplier = Chem.SDMolSupplier('/home/s2451611/MScProject/conformer_dir/5yj0_SMILE.sdf')

for i, mol in enumerate(supplier):
    if mol is not None:
        # Write each conformer to a separate PDB file
        Chem.MolToPDBFile(mol, f'/home/s2451611/MScProject/5yj0/output_{i}.pdb')
