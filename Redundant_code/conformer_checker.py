##############################################################################################
# This script was used to check if each active ligand had at least 10 decoy ligands in order #
# to ensure the data was balanced - redundant as decoys were removed from workflow           #
##############################################################################################

import os
from collections import defaultdict

# Directory path where your files are located
directory_path = '/home/s2451611/MScProject/BROKEN_conformer_dir'

# Initialize a defaultdict to store the count of decoys and actives
file_count = defaultdict(lambda: {'active': 0, 'decoy': 0})

# Iterate over each file in the directory
for file_name in os.listdir(directory_path):
    if file_name.endswith('.sdf'):
        # Split on underscore to get ID and type
        parts = file_name.split('_')
        id = parts[0]
        file_type = 'active' if 'SMILE' in parts[1] else 'decoy'

        # Increment the count in our dictionary
        file_count[id][file_type] += 1

# Check if each ID has 10 times more decoys than actives and write to file if not
with open('missing_decoys.txt', 'w') as f:
    for id, counts in file_count.items():
        if counts['decoy'] < 10 * counts['active']:
            f.write(f'The ID {id} is missing decoys.\n')
