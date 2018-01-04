from alphaspace import *
import mdtraj
import sys
import os


"""
This file loads from a sequence of pdb files created by Dr. Cheng Wang.

The files are in the following format:

il2.apo1.80.33.prot.pdb  {system name}{prot}{.pdb}
il2.apo1.80.33.lig.pdb   {system name}{lig}{.pdb}

"""

def main(files):

    names = set()
    for file in files:
        formatted_name = file.split('.')[:4]
        names.add(".".join(formatted_name))

    print("There are {} snapshots".format(len(names)))

    universe = AS_Universe()

    for name in names:
        universe.set_receptor(structure=mdtraj.load(name+".prot.pdb"),append=True,keepH=True)
        universe.set_binder(structure=mdtraj.load(name+".lig.pdb"),append=True)

    universe.run_alphaspace_mp(10)

    import pickle

    with open("out.as",'wb') as handle:
        pickle.dump(universe,handle)








if __name__ == '__main__':
    main(sys.argv[1:])
