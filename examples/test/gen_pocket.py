"""
This file generate the pocket information and store it in separate dictionary objects.

space.pickl

xyz.pickl

lining.pickl
"""

from alphaspace import *
import mdtraj as md
import sys
import _pickle as pickle




traj = md.load(sys.argv[1])

uni = AS_Universe(traj, guess_receptor_binder=True, guess_by_order=False)

uni.run_alphaspace_mp()

space_dict = {i:{} for i in uni.snapshots_indices}
xyz_dict = {i:{} for i in uni.snapshots_indices}
lining_dict = {i:{} for i in uni.snapshots_indices}

for i in uni.snapshots_indices:
    print("{} of {} generated".format(i,uni.n_frames))


    for j, pocket in enumerate(uni.pockets(i)):

        space_dict[i][j] = pocket.space

        xyz_dict[i][j] = pocket.centroid

        lining_dict[i][j] = pocket.lining_atoms_idx

with open('{}.space.pickl'.format(sys.argv[2]), 'wb') as handle:
    pickle.dump(space_dict, handle)
with open('{}.xyz.pickl'.format(sys.argv[2]), 'wb') as handle:
    pickle.dump(xyz_dict, handle)
with open('{}.lining.pickl'.format(sys.argv[2]), 'wb') as handle:
    pickle.dump(lining_dict, handle)


