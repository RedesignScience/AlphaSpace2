import alphaspace as alpha
import mdtraj as md
import numpy as np
from timeit import default_timer
from collections import defaultdict
import matplotlib.pyplot as plt


def prune_dpockets(d_pocket_dict, sample_ratio=1):
    leader = []
    labels = []
    for d_pocket_idx, pockets in d_pocket_dict.items():
        n_leaders = int(np.ceil(len(pockets) / float(sample_ratio)))
        leader.extend(np.random.choice(pockets, n_leaders, replace=False))
        labels.extend([d_pocket_idx] * n_leaders)
    return leader, labels


traj_path = '/Users/haotian/Desktop/bcl2.holo1.10.h5'
traj = md.load_hdf5(traj_path)
u = alpha.Trajectory(traj, guess_by_order=False, guess_receptor_binder=True)
u.run_alphaspace_mp()

d_pockets = u._gen_d_pockets_iter(sample_frames=20, sample_ratio=1, pocket_space_cutoff=20)
