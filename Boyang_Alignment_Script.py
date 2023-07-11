from rdkit import Chem
from rdkit.Chem import AllChem

# Define the paths to the sdf files
ref_file_path = "/home/s2451611/MScProject/10_percent_crystal_ligand_sdf_db/3u10_ligand.sdf"
conformers_file_path = "/home/s2451611/MScProject/3u10_conf.sdf"

# Load the reference molecule
ref_supplier = Chem.SDMolSupplier(ref_file_path)
reference = next(ref_supplier)

# Load the conformers
conf_supplier = Chem.SDMolSupplier(conformers_file_path)
conformers = [m for m in conf_supplier if m is not None] # Filter out None values

# Align the conformers to the reference
for conformer in conformers:
    AllChem.AlignMol(conformer, reference)

# Save the aligned conformers
writer = Chem.SDWriter('aligned_conformers.sdf')
for conformer in conformers:
    writer.write(conformer)
writer.close()