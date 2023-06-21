from rdkit import Chem
import pandas as pd
import os
import subprocess
from joblib import Parallel, delayed


smile = "CO[C@H]1C[C@H](C)CC2=C(OC)C(=O)C=C(C2=O)NC(=O)/C(=C/C=C[C@@H]([C@H](/C(=C/[C@@H]([C@H]1O)C)/C)OC(=O)N)OC)/C"
output_file = "test.sdf"
command = ["confgen", "-i", smile, "-o", output_file, "--numconf", "20"]
subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
