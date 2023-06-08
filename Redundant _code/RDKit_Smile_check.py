from rdkit import Chem, rdBase, RDLogger
from rdkit.Chem import AllChem

smile = "C[C@H](O)c1cn(N=C234CC[C@]25=C3C(NC(=O)C(C)(C)C)=CC4=C5)cc1NC(=O)N[C@H](C)C1=C(C(=O)O)OC=CC(Cl)=C1"
mol = Chem.MolFromSmiles(smile)
if mol is None:
    print(f"Could not create molecule from SMILES string: {smile}")
else:
    AllChem.EmbedMolecule(mol)
