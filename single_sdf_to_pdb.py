##################################################################################
# This script converts a single sdf file, including all molecules it may contain #
# # and converts each molecule to a pdb file using rdkit                         #
###################################################################################


from rdkit import Chem
from rdkit.Chem import AllChem

# Read the SDF file
supplier = Chem.SDMolSupplier('/home/s2451611/MScProject/Raw_Data/PDBBind_ligand_sdf/1o2o_ligand.sdf', sanitize= False)

for i, mol in enumerate(supplier):
    if mol is not None:
        # Write each conformer to a separate PDB file
        Chem.MolToPDBFile(mol, f'/home/s2451611/MScProject/1yet/output_{i}.pdb')
