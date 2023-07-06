############################################################################################################
# This script copies all pdb files from the PDBbind ligand pdb directory to the active ligand crystal pose #
# pdb bind directory.                                                                                      #
############################################################################################################

import shutil
import os

src_folder = "/home/s2451611/MScProject/Raw_Data/PDBBind_ligand_pdb"
dest_folder = "/home/s2451611/MScProject/Raw_Data/Active_ligand_crystal_poses_pdb_db"

try:
    for filename in os.listdir(src_folder):
        file_path = os.path.join(src_folder, filename)
        if os.path.isfile(file_path):
            shutil.copy(file_path, dest_folder)
    print("All files copied successfully.")
except FileNotFoundError:
    print("Source or destination folder not found.")
except PermissionError:
    print("Permission denied.")
except Exception as e:
    print(f"An error occurred: {e}")
