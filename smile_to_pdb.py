import pandas as pd
import os
from rdkit import Chem, rdBase, RDLogger
from rdkit.Chem import AllChem


def main():
    # Mute warnings
    logger = RDLogger.logger()
    logger.setLevel(RDLogger.CRITICAL)

    # Set the random seed for reproducibility
    rdBase.rdkitRandomSeed = 42

    # Define the input file paths
    active_smiles_path = r"/home/s2451611/MScProject/SCORCH_SMILES/Active_Ligand_Smiles.csv"
    decoy_smiles_path = r"/home/s2451611/MScProject/SCORCH_SMILES/Decoy_Ligand_Smiles.csv"

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
    for index, row in active_df.iterrows():
        pdb_code = row['PDBCode']
        smiles = row['RDKitCanonSmile']
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            AllChem.EmbedMolecule(mol)
            # Create a subdirectory for this PDB code
            pdb_subdir = os.path.join(active_output_dir, pdb_code)
            os.makedirs(pdb_subdir, exist_ok=True)

            # Save the molecule to a PDB file in the subdirectory
            pdb_path = os.path.join(pdb_subdir, f'{pdb_code}_SMILE.pdb')
            AllChem.MolToPDBFile(mol, pdb_path)
        else:
            with open('Broken_Active_Smiles.txt', 'a') as f:
                f.write(f"Invalid SMILES string at index {index}: {smiles}\n")

    print("Processing Decoy SMILES, please wait...")

    # Process the decoy ligands
    for index, row in decoy_df.iterrows():
        ligand = row['Ligand']
        smiles = row['SMILE']
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            AllChem.EmbedMolecule(mol)

            # Save the molecule to a PDB file in the decoy directory
            pdb_path = os.path.join(decoy_output_dir, f'{ligand}_SMILE.pdb')
            AllChem.MolToPDBFile(mol, pdb_path)
        else:
            with open('Broken_Decoy_Smiles.txt', 'a') as f:
                f.write(f"Invalid SMILES string at index {index}: {smiles}\n")

    print("Processing Completed")


if __name__ == "__main__":
    main()
