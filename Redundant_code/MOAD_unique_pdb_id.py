####################################################################################
# This script counts how many unique pdb ids are in the BindingMOAD_2020 directory #
####################################################################################

import os
import re

folder_path = "/home/s2451611/MScProject/Raw_Data/BindingMOAD_2020"

try:
    pdb_ids = set()
    for filename in os.listdir(folder_path):
        match = re.match(r'(\w{4})\..*', filename)
        if match:
            pdb_id = match.group(1)
            pdb_ids.add(pdb_id)
    print(f"There are {len(pdb_ids)} unique PDB IDs.")
except FileNotFoundError:
    print("The directory does not exist.")
except Exception as e:
    print(f"An error occurred: {e}")
