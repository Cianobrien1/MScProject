def convert_pdbqt_to_pdb(pdbqt_file, pdb_file):
    with open(pdbqt_file, 'r') as fin, open(pdb_file, 'w') as fout:
        for line in fin:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                # Write the first 66 characters of the line to the output file
                # This excludes the last column (charge) of the PDBQT file
                fout.write(line[:66] + '\n')
            elif line.startswith('MODEL') or line.startswith('ENDMDL'):
                # Write model delimiters as is
                fout.write(line)

# Example usage:
convert_pdbqt_to_pdb('input.pdbqt', 'output.pdb')
