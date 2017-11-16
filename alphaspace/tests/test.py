import mdtraj

import alphaspace as AS

# initialize universe instance
universe = AS.AS_Universe()
for i in range(1, 1 + 3):
    test_ligand_path = './bcl2/lig/{}.pdb'.format(i)
    test_protein_path = './bcl2/prot/{}.pdb'.format(i)
    ligand = mdtraj.load(test_ligand_path)
    protein = mdtraj.load(test_protein_path)
    universe.set_receptor(protein, append=True)
    universe.set_binder(ligand, append=True)

universe.run_alphaspace_mp()

from scipy.spatial.distance import cdist
import numpy as np

for snapshot_idx in range(universe.n_frames):
    frame_data = universe._data.snapshot(snapshot_idx)

    frame_alpha_idx = frame_data[:, 0].astype(int)
    dist_matrix = cdist(frame_data.xyz(), universe.binder.traj.xyz[snapshot_idx])

    min_idx = np.argmin(dist_matrix, axis=1)
    mins = np.min(dist_matrix, axis=1) * 10  # nm to A
    is_contact = mins < universe.config.hit_dist

    universe._data[frame_alpha_idx, AS.ASDATA_closest_binder_atom_idx] = min_idx
    universe._data[frame_alpha_idx, AS.ASDATA_closest_binder_atom_distance] = mins
    universe._data[frame_alpha_idx, AS.ASDATA_is_contact] = is_contact

for line in universe._data:
    print(line[AS.ASDATA_closest_binder_atom_idx], line[AS.ASDATA_is_contact],
          line[AS.ASDATA_closest_binder_atom_distance])



#
# universe.config.screen_by_lig_cntct = True
# universe.config.screen_by_face = False
#
# universe.screen_pockets()
