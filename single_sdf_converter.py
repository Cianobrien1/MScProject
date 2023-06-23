from rdkit import Chem

# Load a molecule from a .pdb file
mol = Chem.MolFromPDBFile('/home/s2451611/MScProject/Custom_output.pdb')

# Write the molecule to a .sdf file
writer = Chem.SDWriter('Custom_output.sdf')
writer.write(mol)
writer.close()
