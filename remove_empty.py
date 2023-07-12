##################################################
# This script removes empty files in a directory #
##################################################

import os
import argparse

def remove_empty_files(directory):
    for filename in os.listdir(directory):
        file_path = os.path.join(directory, filename)
        # Ensure that the file is actually a file
        if os.path.isfile(file_path):
            # Check if file is empty
            if os.stat(file_path).st_size == 0:
                # Remove the file
                os.remove(file_path)
                print(f"Removed empty file: {file_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Remove empty files in a directory.')
    parser.add_argument('dir', type=str, help='The directory to clean.')
    args = parser.parse_args()
    
    remove_empty_files(args.dir)
