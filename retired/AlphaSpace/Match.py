#AlphaSpace v1.0 - Pocket matching


import os, sys
import AS_fcn


#PDB files to match are entered as command line arguments
pdb_name_list = [arg for arg in sys.argv[1:]]

#input paramters and options are read from the following files
#code will look in cwd first, then in the program directory
param_file = "AS_param.txt"
option_file = "AS_option.txt"


def main():

	AS_fcn.genMatchPocket(pdb_name_list, SAwt = True, metric = 'jaccard',cmethod = 'distance',\
                      isplot = False, isdebug = False)


if __name__ == "__main__":
    main()




