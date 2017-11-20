import sys

from alphaspace import *
import mdtraj


def main(r_path, b_path):
    universe = AS_Universe(receptor=mdtraj.load(r_path), binder=mdtraj.load(b_path))

    """Run alphaspace with 4 cpu cores"""
    universe.run_alphaspace_mp()

    """Iterate over pockets in a particular snapshot"""
    snapshot_idx = 0
    pockets = list(universe.pockets(snapshot_idx))
    pockets.sort(key=lambda p: p.space, reverse=True)

    for pocket in pockets:
        print("Pocket Space {}".format(pocket.space))
        for beta in pocket.betas:
            print("Beta Space {:.2f} at {}".format(beta.space, beta.centroid))


def bcl2_test():
    """Loading 10 snapshots of bcl2 simulation from test/bcl2 folder."""
    universe = AS_Universe()
    for i in range(1, 1 + 10):
        test_ligand_path = './bcl2/lig/{}.pdb'.format(i)
        test_protein_path = './bcl2/prot/{}.pdb'.format(i)
        ligand = mdtraj.load(test_ligand_path)
        protein = mdtraj.load(test_protein_path)
        universe.set_receptor(protein, append=True)
        universe.set_binder(ligand, append=True)
        break

    """Run alphaspace with 4 cpu cores"""
    universe.run_alphaspace_mp()

    for alpha in universe.pockets():


# """Iterate over pockets in a particular snapshot"""
# snapshot_idx = 0
# pockets = list(universe.pockets(snapshot_idx))
# pockets.sort(key=lambda p: p.space, reverse=True)

# for pocket in pockets:
#     for beta in pocket.betas:
#         print(beta.space)



if __name__ == '__main__':
    bcl2_test()
