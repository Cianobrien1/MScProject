from rdkit import Chem
from rdkit.Chem import AllChem

# Read the SDF file
supplier = Chem.SDMolSupplier('/home/s2451611/MScProject/conformer_dir/1yet_SMILE.sdf')

for i, mol in enumerate(supplier):
    if mol is not None:
        # Write each conformer to a separate PDB file
        Chem.MolToPDBFile(mol, f'/home/s2451611/MScProject/1_yet_test/output_{i}.pdb')
