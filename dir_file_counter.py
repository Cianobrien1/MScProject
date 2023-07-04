import os
import sys

"""
This is a Python3 script to count the files and subdirectories in a given directory and its subdirectories.
Use : "python dir_file_counter.py {directory}"
"""

def main():
    def directory_parser(directory):
        file_count = 0
        dir_count = 0

        # Check if the directory exists
        if not os.path.isdir(directory):
            print(f"{directory} is not a valid directory.")
            return

        # Walk through all directories and subdirectories
        for root, dirs, files in os.walk(directory):
            file_count += len(files)
            dir_count += len(dirs)

        print(f"There are {file_count} files in this directory and its subdirectories.")
        print(f"There are {dir_count} subdirectories in this directory.")

    # Check if a directory was provided as a command-line argument
    if len(sys.argv) > 1:
        directory_parser(sys.argv[1])
    else:
        print("Please provide a directory as a command-line argument.")

if __name__ == "__main__":
    main()
