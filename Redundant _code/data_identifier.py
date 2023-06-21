import pandas as pd

def main():
    # Read your CSV file
    df = pd.read_csv('/home/s2451611/MScProject/Summary_by_pdb.csv', low_memory=False)

    # Extract the desired columns
    df = df[["PDB_code", "Source", "Dataset"]]

    # Remove rows beyond 5103
    df = df.iloc[:5103]

    # Save the dataframe to a csv file
    df.to_csv('data_summary.csv', index=False)

    # Remove the 'PDB_' prefix from the 'PDB_code' column
    df['PDB_code'] = df['PDB_code'].apply(lambda x: x.replace('PDB_', '') if isinstance(x, str) else x)

    # Divide the dataframe into subsets based on 'Dataset' and 'Source' and write to txt
    for dataset in df['Dataset'].unique():
        df_dataset = df[df['Dataset'] == dataset]
        df_dataset['PDB_code'].to_csv(f'{dataset}.txt', index=False, header=False)

    for source in df['Source'].unique():
        df_source = df[df['Source'] == source]
        df_source['PDB_code'].to_csv(f'{source}.txt', index=False, header=False)

if __name__ == "__main__":
    main()
