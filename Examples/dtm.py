import mdtraj
from AS_Universe import AS_Universe
from multiprocessing import Pool, cpu_count
import matplotlib
from timeit import default_timer
import numpy as np
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import fcluster, linkage
import matplotlib.pyplot as plt

# initialize universe instance
universe = AS_Universe()
for i in range(1, 1 + 10):
    test_ligand_path = '../Test_system/bcl2/lig/{}.pdb'.format(i)
    test_protein_path = '../Test_system/bcl2/prot/{}.pdb'.format(i)
    ligand = mdtraj.load(test_ligand_path)
    protein = mdtraj.load(test_protein_path)
    universe.set_receptor(protein, append=True)
    universe.set_binder(ligand, append=True)

for i in range(universe.n_frames):
    universe.run(i) # process each frame

universe.config.screen_by_face = True
universe.config.screen_by_ligand_contact = False
universe.screen_pockets()

for d_pocket in universe.d_pockets:
    for polar, nonpolar in d_pocket.pocket_scores:
        print(
            polar, nonpolar
        )
