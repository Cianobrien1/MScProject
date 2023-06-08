import pandas as pd

def missing_pdb_cleaner():
    # Load the data
    decoy_df = pd.read_csv('/home/s2451611/MScProject/SCORCH_SMILES/Decoy_Ligand_Smiles.csv')
    active_df = pd.read_csv('/home/s2451611/MScProject/SCORCH_SMILES/Active_Ligand_Smiles.csv')

    # Extract PDB codes from the decoy dataframe
    decoy_df['PDBCode'] = decoy_df['Ligand'].apply(lambda x: x.split('_')[0])

    # Find the PDB codes that are in the decoy dataframe but not in the active dataframe
    missing_pdb_codes = set(decoy_df['PDBCode']) - set(active_df['PDBCode'])

    # # Write the missing PDB codes to a text file
    # with open('Missing_PDBs.txt', 'w') as f:
    #     for pdb_code in missing_pdb_codes:
    #         f.write(f'{pdb_code}\n')
    
    # Remove rows with missing PDB codes from the decoy dataframe
    decoy_df = decoy_df[~decoy_df['PDBCode'].isin(missing_pdb_codes)]

    # Write the updated decoy dataframe to a new CSV file
    decoy_df.to_csv('/home/s2451611/MScProject/SCORCH_SMILES/Updated_Decoy_Ligand_Smiles.csv', index=False)
    return 

def broken_smiles_cleaner():
    
    """
    Removes the broken SMILES that RDKit could not process stored in the "Broken_Decoy_Smiles.txt" file
    """
    
     # Read the broken ligands
    with open('/home/s2451611/MScProject/SCORCH_SMILES/Template_Broken_Decoy_Smiles.txt', 'r') as f:
        broken_smiles = {line.split()[5].rstrip(":") for line in f}
        
    df = pd.read_csv('/home/s2451611/MScProject/SCORCH_SMILES/Updated_Decoy_Ligand_Smiles.csv')
    
    # Remove rows in the DataFrame that match a broken ligand
    df = df[~df['Ligand'].isin(broken_smiles)]
    
    # Save the filtered DataFrame back to the CSV file
    df.to_csv('/home/s2451611/MScProject/SCORCH_SMILES/Updated_Decoy_Ligand_Smiles.csv', index=False)

def main():
    missing_pdb_cleaner()
    broken_smiles_cleaner()


if __name__ == "__main__":
    main()