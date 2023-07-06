##############################################################################################################
# This script copies the ligands missing from the BindingMOAD dataset, from the pdbqt derrived openbabel pdb #
# directory to a new one containing just the missing ligands                                                 #
##############################################################################################################


import shutil
import os

src_folder = "/home/s2451611/MScProject/openbabel_pdb_crystal_pose_ligands"
dest_folder = "/home/s2451611/MScProject/MOAD_missing_pdbs"
pdb_ids = ["1re8", "1tfz", "1vjy", "2az5", "2pel", "2x2m", "3g5d", "3pvw", "3ri1", "3sfi", "4jzr", "4xrq", 
           "5e95", "5okt", "5uex", "5w7u", "5yj0", "6fp4", "6nlk", "6oa3", "6pgu"]

try:
    # Ensure that the destination directory exists
    if not os.path.exists(dest_folder):
        os.makedirs(dest_folder)

    for pdb_id in pdb_ids:
        file_name = f"{pdb_id}_ligand.pdb"
        src_path = os.path.join(src_folder, file_name)
        dest_path = os.path.join(dest_folder, file_name)

        if os.path.exists(src_path):
            shutil.copy(src_path, dest_path)
        else:
            print(f"File {file_name} does not exist in the source directory.")

    print("Files copied successfully.")
except FileNotFoundError:
    print("Source or destination folder not found.")
except PermissionError:
    print("Permission denied.")
except Exception as e:
    print(f"An error occurred: {e}")
