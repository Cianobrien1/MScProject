##################################################################################################
# This script uses pymol to align the 10 percent test samples to the 10% crystal ligand samples  #
# using sdf file inputs for both ligands being aligned                                           #
##################################################################################################

import os
import glob
from pymol import cmd

# Directories
input_dir1 = "/home/s2451611/MScProject/2023_MMFF_10_percent_sdf"
input_dir2 = "/home/s2451611/MScProject/openbabel_10_percent_sdf_crystal_pose"
output_dir = "MMFF_2023_pymol_aligned_pdbs"

# Create output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

# Output RMSD file
rmsd_file = "MMFF_2023_pymol_avg_rmsd.txt"

def align_and_save():
    all_rmsd_list = []
    for file1 in glob.glob(os.path.join(input_dir1, '*.sdf')):
        pdb_id = os.path.basename(file1).split('_')[0]
        file2 = os.path.join(input_dir2, f"{pdb_id}_ligand.sdf")
        
        if not os.path.exists(file2):
            print(f"No corresponding ligand sdf for {pdb_id}")
            continue
        
        cmd.reinitialize()
        cmd.load(file1, 'conformers')
        cmd.load(file2, 'ligand')
        
        for state in range(1, cmd.count_states('conformers')+1):
            rmsd = cmd.align(f'conformers and state {state}', 'ligand')[0]
            all_rmsd_list.append(rmsd)
        
        # get the state with minimum rmsd
        min_rmsd_state = all_rmsd_list.index(min(all_rmsd_list)) + 1
        cmd.save(os.path.join(output_dir, f"{pdb_id}_aligned.pdb"), f'conformers and state {min_rmsd_state}')
    
    # write average rmsd to file
    with open(rmsd_file, 'w') as f:
        f.write(f"Average RMSD: {sum(all_rmsd_list) / len(all_rmsd_list)}\n")

if __name__ == "__main__":
    align_and_save()
