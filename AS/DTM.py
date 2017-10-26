import mdtraj
from AS import AS_Universe
from multiprocessing import Pool, cpu_count
import matplotlib
from timeit import default_timer
import numpy as np
from scipy.spatial.distance import jaccard,squareform
from scipy.cluster.hierarchy import fcluster,linkage
import matplotlib.pyplot as plt


universe = AS_Universe()
for i in range(1, 1+10):
    test_ligand_path = '../Test_system/bcl2/lig/{}.pdb'.format(1)
    test_protein_path = '../Test_system/bcl2/prot/{}.pdb'.format(1)
    ligand = mdtraj.load(test_ligand_path)
    protein = mdtraj.load(test_protein_path)
    universe.set_receptor(protein, append=True)
    universe.set_binder(ligand, append=True)
    # print(universe.receptor.n_frames)

#
# test_ligand_path = '../Test_system/lig.pdb'
# test_protein_path = '../Test_system/prot_small.pdb'
# universe.set_receptor(mdtraj.load(test_protein_path))
# universe.set_binder(mdtraj.load(test_ligand_path))

universe.config.screen_by_face = False
universe.config.screen_by_ligand_contact = False

for i in range(universe.n_frames):
    universe.run(snapshot_idx=i)

print(sum([universe.cluster(i).n_pockets for i in range(universe.n_frames)]),'total pockets in all frames')


# universe.screen_pockets()
# print(sum([universe.cluster(i).n_pockets for i in range(universe.n_frames)]))
#
# for i in range(universe.n_frames):
#     for pocket in universe.pockets(i):
#         print(pocket.get_lining_atoms())
#     print('end of {}'.format(i))

pockets = []
for i in range(universe.n_frames):
    pockets.extend(list(universe.cluster(snapshot_idx=i).pockets))



pocket_diff = np.zeros((len(pockets),len(pockets)))

lining_atom_set = set()
for pocket in pockets:
    lining_atom_set = lining_atom_set.union(pocket.get_lining_atoms())

for i in range(len(pockets)):
    for j in range(i,len(pockets)):
        p1 = set(pockets[i].get_lining_atoms())
        p2 = set(pockets[j].get_lining_atoms())
        diff = p1.symmetric_difference(p2)
        union = p1.union(p2)
        pocket_diff[i,j] = pocket_diff[j,i] = len(diff)/len(union)


pocket_diff = squareform(pocket_diff)
_linkage = linkage(pocket_diff,method='complete')
clustered_list = list(fcluster(Z=_linkage,t=0.75,criterion='distance'))


# print(pocket_diff.shape)
#
# plt.imshow(pocket_diff, cmap=plt.cm.gray, interpolation='nearest')
# plt.show()






