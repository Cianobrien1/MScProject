########################################################################################
# This script was used to check how many .bio1 files were in the BindingMOAD directory #
########################################################################################

import os

count = 0

for file in os.listdir("/home/s2451611/MScProject/BindingMOAD_2020"):
    if file.endswith(".bio1"):
        count +=1 
    else: 
        continue

print(count)