###################################################################################
# This script was used to convert all pdbqt files to pdb files in a custom manner #
# it was replaced by using openbabel due to errors in the converion               #
###################################################################################

from joblib import Parallel, delayed
import os
import glob

def convert_pdbqt_to_pdb(pdbqt_file, pdb_file):
    print(f"Converting {pdbqt_file} to {pdb_file}...")
    with open(pdbqt_file, 'r') as fin, open(pdb_file, 'w') as fout:
        for line in fin:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                fout.write(line[:66] + '\n')
            elif line.startswith('MODEL') or line.startswith('ENDMDL'):
                fout.write(line)
    print(f"Conversion of {pdbqt_file} completed.")

input_dir = "/home/s2451611/MScProject/all_receptors_ligands/all_receptors_ligands"
output_dir = "/home/s2451611/MScProject/pdb_all_receptor_ligands"

# Get a list of all PDBQT files in the input directory and its subdirectories
pdbqt_files = glob.glob(os.path.join(input_dir, "**/*.pdbqt"), recursive=True)

print(f"Found {len(pdbqt_files)} PDBQT files.")

# Use joblib to convert the files in parallel
for pdbqt_file in pdbqt_files:
    # Construct the output file path
    relative_path = os.path.relpath(pdbqt_file, input_dir)
    pdb_file = os.path.join(output_dir, relative_path.replace('.pdbqt', '.pdb'))

    # Create the output directory if it doesn't exist
    os.makedirs(os.path.dirname(pdb_file), exist_ok=True)

    # Convert the file
    convert_pdbqt_to_pdb(pdbqt_file, pdb_file)
