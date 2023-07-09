
########################################################################
# This script aligns a single pdb file to another pdb file using rkdit #
########################################################################
from rdkit import Chem
from rdkit.Chem import rdMolAlign

def align_pdb_files(ref_pdb_file, fit_pdb_file, output_pdb_file):
    # Load the reference and fit molecules
    ref_mol = Chem.MolFromPDBFile(ref_pdb_file)
    fit_mol = Chem.MolFromPDBFile(fit_pdb_file)

    # Perform the alignment
    rmsd = rdMolAlign.AlignMol(fit_mol, ref_mol)

    # Print the RMSD
    print('RMSD:', rmsd)

    # Save the aligned molecule
    writer = Chem.PDBWriter(output_pdb_file)
    writer.write(fit_mol)
    writer.close()

# Specify the reference PDB file, fit PDB file, and output PDB file
ref_pdb_file = '/home/s2451611/MScProject/10_percent_crystal_ligand_pdb_db/3u10_ligand.pdb'
fit_pdb_file = '/home/s2451611/MScProject/pdb_10_percent_alignment_test_sample/3u10_SMILE_18.pdb'
output_pdb_file = 'aligned.pdb'

align_pdb_files(ref_pdb_file, fit_pdb_file, output_pdb_file)
