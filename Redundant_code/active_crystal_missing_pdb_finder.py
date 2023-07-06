#####################################################################################################
 # This script was used to find out what 2 pdb id codes were missing from the Active ligand crytsal #
 # pose directory, (database derrived)                                                              #
#####################################################################################################

import os
import re

folder1 = "/home/s2451611/MScProject/conformer_dir"
folder2 = "/home/s2451611/MScProject/Raw_Data/Active_ligand_crystal_poses_pdb_db"

try:
    pdb_ids1 = set()
    for filename in os.listdir(folder1):
        match = re.match(r'(\w{4})_SMILE.sdf', filename)
        if match:
            pdb_ids1.add(match.group(1))

    pdb_ids2 = set()
    for filename in os.listdir(folder2):
        match = re.match(r'(\w{4})_ligand.pdb', filename)
        if match:
            pdb_ids2.add(match.group(1))

    missing_ids = pdb_ids1 - pdb_ids2
    print(f"The missing PDB IDs are: {missing_ids}")
except FileNotFoundError:
    print("One of the directories does not exist.")
except Exception as e:
    print(f"An error occurred: {e}")
