##################################################################################
# This was an experimental way to fix issues with SMILES I was having by using a #
# SMILE cleaner from vincent's github"                                           #
##################################################################################

from rdkit import Chem
import re

def smiles_cleaner(smiles, return_idx=False, deep_clean=True):
    smiles = list(smiles)
    idx = []
    idx_bad = []
    clean_smiles = []
    for i, s in enumerate(smiles):
        try:
            if deep_clean:
                smi = s.replace(" ", "")
                smi = max(smi.split("."), key=len)
                for m in re.findall("\[.*?\]", smi):
                    if m not in ['[N+]', '[O-]']:
                        smi = smi.replace(m, m[1].upper())
                smi = smi.replace("@", "")
                smi = smi.replace("/C", "C")
                smi = smi.replace("\\C", "C")
                smi = smi.replace("/c", "c")
                smi = smi.replace("\\c", "c")
                
            else:
                smi = s
            
            m = Chem.MolFromSmiles(smi, sanitize=True)
        except:
            m = None

        if m is not None:
            idx.append(i)
            clean_smiles.append(smi)
        else:
            idx_bad.append(i)
            print(f"Warning: invalid SMILES in position {i}: {s}")

    if return_idx:
        return clean_smiles, idx, idx_bad
    else:
        return clean_smiles

# The SMILES string you want to clean
smiles_string = "CSSSS[S@H]1S23(F)(=S145(SF)[C@H](S)CS[S@@]24[S@H]5SF)S[S@H]3C=S"

# Call the function with the SMILES string
cleaned_smiles = smiles_cleaner([smiles_string], deep_clean=True)

# Print the cleaned SMILES string
print(cleaned_smiles)
