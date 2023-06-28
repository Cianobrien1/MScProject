###################################################################
# This script is used to count the number of files in a directory # 
###################################################################

import sys
import os

"""
This is a simple python3 script to count the files in a given directory
Use : "python dir_file_counter.py {directory}"

"""

def main():
    def file_counter(directory):
        # Check if the directory exists
        if not os.path.isdir(directory):
            print(f"{directory} is not a valid directory.")
            return

        # Use a generator to count the files
        count = sum(1 for entry in os.scandir(directory) if entry.is_file())

        print(f"There are {count} files in this directory.")
        return 

    # Check if a directory was provided as a command-line argument
    if len(sys.argv) > 1:
        file_counter(sys.argv[1])
    else:
        print("Please provide a directory as a command-line argument.")

if __name__ == "__main__":
    main()
