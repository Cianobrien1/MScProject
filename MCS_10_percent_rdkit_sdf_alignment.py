##################################################################################################
# This script uses an MCS approach and rdkit to align the 10 percent test samples                #
# to the 10% crystal ligand samples using sdf file inputs for both ligands being aligned         #
##################################################################################################

import os
from rdkit import Chem, RDConfig
from rdkit.Chem import AllChem, rdFMCS
from joblib import Parallel, delayed
from tqdm import tqdm

crystal_dir = "/home/s2451611/MScProject/openbabel_10_percent_sdf_crystal_pose"
conformers_dir = "/home/s2451611/MScProject/SANITIZED_10_percent_conformer_dir"
output_dir = "/home/s2451611/MScProject/MCS_openbabel_aligned_pdbs"
rmsd_file = "/home/s2451611/MScProject/MCS_openbabel_rdkit_conformer_rmsd.txt"
fail_log_file = "/home/s2451611/MScProject/failed_alignment.txt"


def process_pdb_id(pdb_id):
    # Read the crystal pose and the conformers
    crystal_pose_file = os.path.join(crystal_dir, f"{pdb_id}_ligand.sdf")
    crystal_pose_supplier = Chem.SDMolSupplier(crystal_pose_file)
    if not crystal_pose_supplier:
        return []
    crystal_pose = crystal_pose_supplier[0]
    if crystal_pose is None:
        return []

    conformers_file = os.path.join(conformers_dir, f"{pdb_id}_SMILE.sdf")
    conformers_supplier = Chem.SDMolSupplier(conformers_file)
    if not conformers_supplier:
        return []
    conformers = [m for m in conformers_supplier if m is not None]
    if not conformers:
        return []

    # Find the MCS
    mcs = rdFMCS.FindMCS([crystal_pose] + list(conformers), 
                         threshold=0.8, 
                         completeRingsOnly=True, 
                         ringMatchesRingOnly=True)

    # Align each conformer to the crystal pose
    patt = Chem.MolFromSmarts(mcs.smartsString)
    refMol = crystal_pose
    refMatch = refMol.GetSubstructMatch(patt)
    
    rmsVs = []
    for probeMol in conformers:
        try:
            mv = probeMol.GetSubstructMatch(patt)
            rms = AllChem.AlignMol(probeMol, refMol, atomMap=list(zip(mv, refMatch)))
            rmsVs.append(rms)
        except RuntimeError:
            with open(fail_log_file, 'a') as fail_log:
                fail_log.write(f"No common substructure found between the probe and reference molecule for PDB ID {pdb_id}\n")

    if rmsVs: # Check if the list is not empty
        # Save the conformer with the lowest RMSD to a pdb file
        min_rmsd_index = rmsVs.index(min(rmsVs))  # Get the index of the conformer with lowest RMSD
        best_conformer = conformers[min_rmsd_index]
        
        writer = Chem.PDBWriter(os.path.join(output_dir, f"{pdb_id}_best.pdb"))
        writer.write(best_conformer)
        writer.close()

    return rmsVs

pdb_ids = [filename.split('_')[0] for filename in os.listdir(crystal_dir)]

results = Parallel(n_jobs=-1)(delayed(process_pdb_id)(pdb_id) for pdb_id in tqdm(pdb_ids))

all_rmsds = [rms for sublist in results for rms in sublist]

# Calculate average RMSD across all conformers and save to file
avg_rmsd = sum(all_rmsds) / len(all_rmsds)
with open(rmsd_file, 'w') as f:
    f.write(f"Average RMSD: {avg_rmsd}\n")
