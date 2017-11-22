import sys

from alphaspace import *
import mdtraj


def main(top, coord):
    """Loading 10 snapshots of bcl2 simulation from test/bcl2 folder."""
    traj = mdtraj.load(coord, top=top)

    traj.remove_solvent(inplace=True)

    non_h_atoms = []
    for atom in traj.top.atoms:
        if atom.element.number != 0:
            non_h_atoms.append(atom.index)
    traj.atom_slice(np.array(non_h_atoms), inplace=True)

    receptor_idx = []
    binder_idx = []

    for residue in traj.top.residues:
        if residue.index < 172:
            receptor_idx.extend([atom.index for atom in residue.atoms])
        elif residue.index < 204:
            binder_idx.extend([atom.index for atom in residue.atoms])

    receptor_traj = traj.atom_slice(np.array(receptor_idx))
    binder_traj = traj.atom_slice(np.array(binder_idx))

    universe = AS_Universe(receptor=receptor_traj, binder=binder_traj)

    universe.config.screen_by_lig_cntct = True
    """Run alphaspace with 4 cpu cores"""
    universe.run_alphaspace_mp()

    for frame in range(universe.n_frames):
        for atom in universe.binder.atoms:
            atom.linked_alpha = []

        for pocket in universe.pockets(frame, active_only=True):
            for alpha in pocket.alphas:
                universe.binder.atom(alpha.closest_atom_idx).linked_alpha.append(alpha)

        total_space = 0.0
        occupied_space = 0.0
        for atom in universe.binder.atoms:
            if (len(atom.linked_alpha)) > 0:
                total_space += sum([alpha.space for alpha in atom.linked_alpha])
                occupied_space += sum([alpha.space for alpha in atom.linked_alpha if alpha.is_contact])

        print(total_space, occupied_space)


# """Iterate over pockets in a particular snapshot"""
# snapshot_idx = 0
# pockets = list(universe.pockets(snapshot_idx))
# pockets.sort(key=lambda p: p.space, reverse=True)

# for pocket in pockets:
#     for beta in pocket.betas:
#         print(beta.space)



if __name__ == '__main__':
    coord = sys.argv[1]
    top = sys.argv[2]

    main(top, coord)
