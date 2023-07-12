#############################################################################
# This script is used to generata RDKit 2023 SMILES for all pdb files input #
#############################################################################

from rdkit import Chem
from rdkit.Chem import AllChem
import pandas as pd
import glob
import os

def pdb_to_smiles(file_path):
    try:
        # Load the PDB file as a RDKit Mol object
        mol = Chem.MolFromPDBFile(file_path, removeHs=False)

        # Generate the canonical SMILES string
        smiles = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
    except Exception as e:
        print(f"Error converting {file_path} to SMILES: {e}")
        smiles = None

    # Extract the PDB ID from the file name
    pdb_id = os.path.basename(file_path).split('_')[0]

    return pdb_id, smiles

def main():
    # Define the directory containing the PDB files
    pdb_dir = "/home/s2451611/MScProject/openbabel_pdb_crystal_pose_ligands"

    # Define the output CSV file path
    output_file_path = "/home/s2451611/MScProject/SCORCH_SMILES/RDKIT_2023_Active_SMILES.csv"

    # Find all PDB files in the directory
    pdb_files = glob.glob(os.path.join(pdb_dir, "*.pdb"))

    # Convert each PDB file to a SMILES string
    data = [pdb_to_smiles(file_path) for file_path in pdb_files]

    # Create a pandas DataFrame from the data
    df = pd.DataFrame(data, columns=['PDBCode', 'RDKitCanonSmile'])

    # Save the DataFrame to a CSV file
    df.to_csv(output_file_path, index=False)

if __name__ == "__main__":
    main()
