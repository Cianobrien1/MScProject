import pandas as pd

def main():
    # Load the data
    decoy_df = pd.read_csv('/home/s2451611/MScProject/SCORCH_SMILES/Decoy_Ligand_Smiles.csv')
    active_df = pd.read_csv('/home/s2451611/MScProject/SCORCH_SMILES/Active_Ligand_Smiles.csv')

    # Extract PDB codes from the decoy dataframe
    decoy_df['PDBCode'] = decoy_df['Ligand'].apply(lambda x: x.split('_')[0])

    # Find the PDB codes that are in the decoy dataframe but not in the active dataframe
    missing_pdb_codes = set(decoy_df['PDBCode']) - set(active_df['PDBCode'])

    # Write the missing PDB codes to a text file
    with open('Missing_PDBs.txt', 'w') as f:
        for pdb_code in missing_pdb_codes:
            f.write(f'{pdb_code}\n')


if __name__ == "__main__":
    main()