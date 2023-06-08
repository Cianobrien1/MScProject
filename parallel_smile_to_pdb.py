import pandas as pd
import os
from rdkit import Chem, rdBase, RDLogger
from rdkit.Chem import AllChem
from joblib import Parallel, delayed
    
def process_smiles(row, output_dir, is_decoy=False):
    # Mute warnings
    rdBase.DisableLog('rdApp.*')
    
    # Process the ligands
    if is_decoy:
        ligand = row['Ligand']
        smiles = row['SMILE']
        filename = 'Broken_Decoy_Smiles.txt'
    else:
        ligand = row['PDBCode']
        smiles = row['RDKitCanonSmile']
        filename = 'Broken_Active_Smiles.txt'

    mol = Chem.MolFromSmiles(smiles)

    if mol is not None:
        # Add temporary hydrogens & remove before saving to pdb file
        h_mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(h_mol)
        mol = Chem.RemoveHs(h_mol)

        # Save the molecule to a PDB file in the directory
        pdb_path = os.path.join(output_dir, f'{ligand}_SMILE.pdb')
        AllChem.MolToPDBFile(mol, pdb_path)
    else:
        with open(filename, 'a') as f:
            f.write(f"Invalid SMILES string for ligand {ligand}: {smiles}\n")


def main():

    # Set the random seed for reproducibility
    rdBase.rdkitRandomSeed = 42

    # Define the input file paths
    active_smiles_path = r"/home/s2451611/MScProject/SCORCH_SMILES/Active_Ligand_Smiles.csv"
    decoy_smiles_path = r"/home/s2451611/MScProject/SCORCH_SMILES/Updated_Decoy_Ligand_Smiles.csv"

    # Define the output directory paths
    active_output_dir = r"/home/s2451611/MScProject/SMILE_pdbs"
    decoy_output_dir = r"/home/s2451611/MScProject/Decoy_Smile_PDBs"

    # Ensure the output directories exist
    os.makedirs(active_output_dir, exist_ok=True)
    os.makedirs(decoy_output_dir, exist_ok=True)

    # Read the CSV files
    active_df = pd.read_csv(active_smiles_path, dtype={'PDBCode': str})
    decoy_df = pd.read_csv(decoy_smiles_path, dtype={'Ligand': str})

    print("Processing Active SMILES, please wait...")

    # Process the active ligands
    Parallel(n_jobs=-1)(delayed(process_smiles)(row, active_output_dir, False) for _, row in active_df.iterrows())

    print("Processing Decoy SMILES, please wait...")

    # Process the decoy ligands
    Parallel(n_jobs=-1)(delayed(process_smiles)(row, decoy_output_dir, True) for _, row in decoy_df.iterrows())

    print("Processing Completed")


if __name__ == "__main__":
    main()
