from rdkit import Chem
import os
from joblib import Parallel, delayed

output_dir = "/home/s2451611/MScProject/conformer_dir/"
input_dir = "/home/s2451611/MScProject/pdb_to_mols_dir/"

def confgen(file_arg):
    output_file_name = os.path.basename(file_arg).rstrip(".mol") + ".sdf"
    output_file = output_dir + output_file_name
    os.system(f"confgen -i {file_arg} -o {output_file}")

def main():
    os.makedirs(output_dir, exist_ok=True)
    file_names = os.listdir(input_dir)
    Parallel(n_jobs=-1)(delayed(confgen)(os.path.join(input_dir, file_name)) for file_name in file_names)

if __name__ == "__main__":
    main()
