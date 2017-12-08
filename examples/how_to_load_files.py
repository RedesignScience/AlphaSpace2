"""
Below is an example of how to load structure files into alphaspace universe container.
You can either initialize AS_Universe instance with trajectory object,
or add in the trajectory later with the set_binder,set_receptor method
"""

import alphaspace
import mdtraj

"""load 1 snapshot in PDB format"""

receptor = mdtraj.load('./bcl2/prot/1.pdb')
binder = mdtraj.load('./bcl2/lig/1.pdb')
universe = alphaspace.AS_Universe(receptor=receptor, binder=binder)

universe.run_alphaspace_mp()

"""Load 10 snapshot in PDB format"""
universe = alphaspace.AS_Universe()
for i in range(1, 1 + 10):
    test_ligand_path = './bcl2/lig/{}.pdb'.format(i)
    test_protein_path = './bcl2/prot/{}.pdb'.format(i)
    ligand = mdtraj.load(test_ligand_path)
    protein = mdtraj.load(test_protein_path)
    universe.set_receptor(protein, append=True)
    universe.set_binder(ligand, append=True)

universe.run_alphaspace_mp()

"""Load several snapshot in a trajectory file with topology file"""
trajectory = mdtraj.load("Path to Trajectory", top="Path to Topology")

universe = alphaspace.AS_Universe(receptor=trajectory, guess_receptor_binder=True, guess_by_order=True)

universe.run_alphaspace()
