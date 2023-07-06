##########################################################################
# This script counts the lines in given input file and prints the number #
##########################################################################

import argparse

def count_lines_efficient(file_path):
    try:
        with open(file_path, 'r') as file:
            count = sum(1 for _ in file)
            return count
    except FileNotFoundError:
        return "The file does not exist."

parser = argparse.ArgumentParser(description="Count the number of lines in a file.")
parser.add_argument("--file", "-f", help="The file to count lines in.")
args = parser.parse_args()

print(f"The file has {count_lines_efficient(args.file)} lines.")
