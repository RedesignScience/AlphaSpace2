from alphaspace import *
import mdtraj
import sys
import os
import numpy as np

import _pickle as pickle

"""
This file loads from a sequence of pdb files created by Dr. Cheng Wang.

The files are in the following format:

il2.apo1.80.33.prot.pdb  {system name}{prot}{.pdb}
il2.apo1.80.33.lig.pdb   {system name}{lig}{.pdb}

"""


def run_alphaspace(files):
    """Holo naming  il2.holo1.80.76.protlig.pdb"""
    traj = mdtraj.load(files)

    protein_atoms = np.array([atom.index for atom in traj.top.atoms if atom.residue.index < 128])
    ligand_atoms = np.array([atom.index for atom in traj.top.atoms if atom.residue.index >= 128])
    protein = traj.atom_slice(np.array(protein_atoms))
    ligand = traj.atom_slice(np.array(ligand_atoms))

    universe = AS_Universe(receptor=protein, binder=ligand)

    print(universe)

    """Apo Naming Scheme"""
    # names = set()
    # for file in files:
    #     formatted_name = file.split('.')[:4]
    #     names.add(".".join(formatted_name))
    #
    # print("There are {} snapshots".format(len(names)))
    #
    # universe = AS_Universe()
    # for name in names:
    #     universe.set_receptor(structure=mdtraj.load(name + ".prot.pdb"), append=True, keepH=True)
    #     universe.set_binder(structure=mdtraj.load(name + ".lig.pdb"), append=True)

    universe.run_alphaspace_mp(4)


    with open("out.as", 'wb') as handle:
        pickle.dump(universe, handle)


def write_d_pockets(as_file):


    with open(as_file,'rb') as handle:
        universe = pickle.load(handle)

    universe.config = AS_Config()


    print("fluctuation  space   space_std  occupancy occurrence   total_snapshot")
    for d_pocket in universe.d_pockets:
        xyz = np.array([pocket.centroid for pocket in d_pocket])
        centroid = np.average(xyz,axis=0)
        fluctuation = np.average(
                                [np.linalg.norm(coord - centroid) for coord in xyz])
        spaces = np.array([pocket.space for pocket in d_pocket])
        space = np.average(spaces)
        space_std = np.std(spaces)

        occurrence = len(set([pocket.snapshot_idx for pocket in d_pocket]))

        occupancy = np.average([pocket.occupancy for pocket in d_pocket])

        total_snapshot = universe.n_frames

        print(np.round(fluctuation,2),np.round(space, 1) ,np.round(space_std,1) ,np.round(occupancy,3),occurrence, total_snapshot)








if __name__ == '__main__':
    # run_alphaspace(sys.argv[1:])
    write_d_pockets(sys.argv[1])
