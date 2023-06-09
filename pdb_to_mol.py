from rdkit import Chem
import os
from joblib import Parallel, delayed

output_dir = "/home/s2451611/MScProject/pdb_to_mols_dir"
os.makedirs(output_dir, exist_ok=True)

def pdb_to_mol(pdb_file_path):
    molecule = Chem.rdmolfiles.MolFromPDBFile(pdb_file_path, sanitize=False)
    
    # # These lines can be activated if Chem.rdmolfiles.MolFromPDBFile() is returning a None type to find which pdb files are causing an issue
    # if molecule is None:
    #     with open("broken_pdb.txt", "a")as f:
    #         f.write(f"Could not convert {pdb_file_path} to a molecule.\n")
    # else:
    mol_file = os.path.basename(pdb_file_path).rstrip(".pdb")
    Chem.rdmolfiles.MolToMolFile(molecule, f"{output_dir}/{mol_file}.mol")

def main():
    input_dirs = ["/home/s2451611/MScProject/Acitve_Smiles_pdbs", "/home/s2451611/MScProject/Decoy_Smiles_pdbs"]

    for input_dir in input_dirs:
        pdb_files = [os.path.join(input_dir, f) for f in os.listdir(input_dir)]
        Parallel(n_jobs=-1)(delayed(pdb_to_mol)(pdb_file) for pdb_file in pdb_files)

if __name__ == "__main__":
    main()




