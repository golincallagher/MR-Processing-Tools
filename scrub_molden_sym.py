"""
scrub_molden_sym.py
--------------------
Removes 'Sym=' lines from a MOLPRO-generated .molden file.

MOLPRO writes a 'Sym=' label on each molecular orbital line in .molden
output. Some visualization software (e.g. VMD) cannot parse .molden
files containing this label, causing orbital data to be unreadable.
This script strips all 'Sym=' lines and writes the result to a new file,
leaving the original file untouched.

Usage:
    python scrub_molden_sym.py <input.molden> <output.molden>

Example:
    python scrub_molden_sym.py edgn2n2coms5cas1.molden edgn2n2coms5casorbitals1.molden
"""

import sys


def remove_sym_lines(input_file, output_file):
    """
    Reads input_file, removes every line containing 'Sym=', and writes
    the remaining lines to output_file.
    """
    with open(input_file, 'r') as f:
        lines = f.readlines()

    filtered_lines = [line for line in lines if 'Sym=' not in line]

    with open(output_file, 'w') as f:
        f.writelines(filtered_lines)


def main():
    if len(sys.argv) != 3:
        print("Usage: python scrub_molden_sym.py <input.molden> <output.molden>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    remove_sym_lines(input_file, output_file)
    print(f"Done. Wrote cleaned file to '{output_file}'.")


if __name__ == '__main__':
    main()
