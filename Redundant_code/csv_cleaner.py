#############################################################################
# This script was used to "clean" the active and decoy ligand csv files by  #
# removing broken SMILES and any active with a decoy that had broken SMILES #
#############################################################################

import pandas as pd

def missing_pdb_cleaner():
    # Load the data
    decoy_df = pd.read_csv('/home/s2451611/MScProject/SCORCH_SMILES/Decoy_Ligand_Smiles.csv')
    active_df = pd.read_csv('/home/s2451611/MScProject/SCORCH_SMILES/Active_Ligand_Smiles.csv')

    # Extract PDB codes from the decoy dataframe
    decoy_df['PDBCode'] = decoy_df['Ligand'].apply(lambda x: x.split('_')[0])

    # Find the PDB codes that are in the decoy dataframe but not in the active dataframe
    missing_pdb_codes = set(decoy_df['PDBCode']) - set(active_df['PDBCode'])
    
    # Remove rows with missing PDB codes from the decoy dataframe
    decoy_df = decoy_df[~decoy_df['PDBCode'].isin(missing_pdb_codes)]

    # Write the updated decoy dataframe to a new CSV file
    decoy_df.to_csv('/home/s2451611/MScProject/SCORCH_SMILES/Updated_Decoy_Ligand_Smiles.csv', index=False)
    return 

def broken_smiles_and_active_cleaner():
    """
    Removes the broken SMILES that RDKit could not process stored in the "Broken_Decoy_Smiles.txt" file
    and also updates the active ligands dataframe by removing rows with matching PDB codes.
    """
    
    # Read the broken ligands
    with open('/home/s2451611/MScProject/SCORCH_SMILES/Template_Broken_Decoy_Smiles.txt', 'r') as f:
        broken_pdb_codes = {line.split()[5].split('_')[0].rstrip(":") for line in f}

    df_decoy = pd.read_csv('/home/s2451611/MScProject/SCORCH_SMILES/Updated_Decoy_Ligand_Smiles.csv')
    
    # Remove rows in the Decoy DataFrame that match a broken ligand
    df_decoy = df_decoy[~df_decoy['Ligand'].isin(broken_pdb_codes)]
    
    # Save the filtered DataFrame back to the Decoy CSV file
    df_decoy.to_csv('/home/s2451611/MScProject/SCORCH_SMILES/Updated_Decoy_Ligand_Smiles.csv', index=False)

    # Load the active ligands dataframe
    df_active = pd.read_csv('/home/s2451611/MScProject/SCORCH_SMILES/Active_Ligand_Smiles.csv')
    
    # Remove rows in the Active DataFrame that match a broken PDB code
    df_active = df_active[~df_active['PDBCode'].isin(broken_pdb_codes)]
    
    # Save the filtered DataFrame back to the Active CSV file
    df_active.to_csv('/home/s2451611/MScProject/SCORCH_SMILES/Updated_Active_Ligand_Smiles.csv', index=False)
    
    df_decoy = pd.read_csv('/home/s2451611/MScProject/SCORCH_SMILES/Updated_Decoy_Ligand_Smiles.csv')
    
    # Remove rows in the Decoy DataFrame that match a broken PDB code
    df_decoy = df_decoy[~df_decoy['PDBCode'].isin(broken_pdb_codes)]
    
    # Save the filtered DataFrame back to the Decoy CSV file
    df_decoy.to_csv('/home/s2451611/MScProject/SCORCH_SMILES/Updated_Decoy_Ligand_Smiles.csv', index=False)

def main():
    missing_pdb_cleaner()
    broken_smiles_and_active_cleaner()


if __name__ == "__main__":
    main()
