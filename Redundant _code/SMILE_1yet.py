####################################################################################
# This script was used to generate a SMILE from thr 1yet mol2 file as the original #
# SMILE provided for this ligand could not be read by RDKit                        #
####################################################################################


from rdkit import Chem

input_file = "/home/s2451611/MScProject/1yet/1yet_ligand.mol2"

def generate_smiles_from_mol2(filename):
    suppl = Chem.MolFromMol2File(filename, sanitize=False)
    smiles = Chem.MolToSmiles(suppl)
    return smiles

def save_to_file(smile):
    with open("output.smi", "w") as f:
        f.write(f"{smile}\n")
    return

def main():
    save_to_file(generate_smiles_from_mol2(input_file))
    
if __name__ == "__main__":
    main()