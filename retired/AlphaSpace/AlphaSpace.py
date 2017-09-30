#!~/anaconda3/bin/python
#AlphaSpace v1.0 - FCTM


import sys

from retired.AlphaSpace import AS_fcn

#PDB file(s) to map is/are entered as command line argument (1 file if protein and ligand together...
# 2 files if protein (first) and ligand (second) separate)
pdb_file_name = [f_name for f_name in sys.argv[1:]]

if len(sys.argv[1:]) > 1:
	pdb_lines = AS_fcn.combine_prot_lig(sys.argv[1:])
else:
	pdb_lines = [line for line in open(sys.argv[1])]

#input paramters and options are read from the following files
#code will look in cwd first, then in the program directory
param_file = "AS_param.txt"
option_file = "AS_option.txt"


def main():

	pdb, ss = AS_fcn.runAlphaSpace(pdb_lines, default = False, isWrite=True)


if __name__ == "__main__":
    main()




