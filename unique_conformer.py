import os


conformer_dir = "/home/s2451611/MScProject/BROKEN_conformer_dir"

# List comprehension to get IDs from file names
# os.listdir(directory_path) returns a list of all file names in the directory
# file_name.split('_')[0] splits each file name on the underscore and takes the first part (the ID)
# The set(...) part removes duplicates
ids = set(file_name.split('_')[0] for file_name in os.listdir(conformer_dir) if file_name.endswith('.sdf'))

# Now, we'll write these IDs to a text file
with open('conformer_pdbcodes.txt', 'w') as f:
    for id in ids:
        f.write(id + '\n')
