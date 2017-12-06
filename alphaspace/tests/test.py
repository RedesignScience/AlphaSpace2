import sys

from alphaspace import *
import mdtraj
import networkx
from itertools import combinations, chain
from AS_Funct import combination_intersection_count, combination_union_count


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

    universe.config.screen_by_lig_cntct = True
    """Run alphaspace with 4 cpu cores"""
    universe.run_alphaspace_mp()

    for atom in universe.binder.atoms:
        atom.linked_alpha = []

    for pocket in universe.pockets(active_only=True):
        for alpha in pocket.alphas:
            universe.binder.atom(alpha.closest_atom_idx).linked_alpha.append(alpha)

    for atom in universe.binder.atoms:
        if (len(atom.linked_alpha)) > 0:
            print(len(atom.linked_alpha))


def generate_communities():
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

    # universe.config.screen_by_lig_cntct = True
    """Run alphaspace with 4 cpu cores"""
    universe.run_alphaspace_mp()

    universe._gen_communities()

    for community in universe.communities[0]:
        print('community')
        print(len(community))
        print(len([p for p in community if p.core_aux_minor == 'core']))
        print(len([p for p in community if p.core_aux_minor == 'aux']))
        print(len([p for p in community if p.core_aux_minor == 'minor']))


def generate_dtm_communities():
    universe = AS_Universe()
    for i in range(1, 1 + 10):
        test_ligand_path = './bcl2/lig/{}.pdb'.format(i)
        test_protein_path = './bcl2/prot/{}.pdb'.format(i)
        ligand = mdtraj.load(test_ligand_path)
        protein = mdtraj.load(test_protein_path)
        universe.set_receptor(protein, append=True)
        universe.set_binder(ligand, append=True)

    universe.config.screen_by_lig_cntct = True
    universe.config.screen_by_face = False

    universe.run_alphaspace_mp()
    universe._gen_communities()

    #  find all core d pockets
    d_pockets = universe._gen_d_pockets()
    core_d_pockets = {}
    for i, pockets in d_pockets.items():
        for pocket in pockets:
            if pocket.core_aux_minor == "core":
                core_d_pockets[i] = pockets
                for p in pockets:
                    p._is_core_d = True
                break

    # mark all pockets connected to members of d-pocket
    for snapshot_idx in universe.snapshots_indices:
        network = universe._pocket_network[snapshot_idx]
        for p in network.nodes():
            if p._is_core_d:
                for cp in network[p]:
                    cp._connected = True

    for i, pockets in d_pockets.items():
        if i not in core_d_pockets:
            print(len([True for p in pockets if p._connected]))










if __name__ == '__main__':
    generate_dtm_communities()
