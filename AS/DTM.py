import mdtraj
from AS import AS_Universe
from multiprocessing import Pool, cpu_count
import matplotlib
from timeit import default_timer

universe = AS_Universe()
for i in range(1, 1+3):
    test_ligand_path = '../Test_system/bcl2/lig/{}.pdb'.format(i)
    test_protein_path = '../Test_system/bcl2/prot/{}.pdb'.format(i)
    ligand = mdtraj.load(test_ligand_path)
    protein = mdtraj.load(test_protein_path)
    universe.set_receptor(protein, append=True)
    universe.set_binder(ligand, append=True)
    # print(universe.receptor.n_frames)

universe.config.screen_by_face = True
universe.config.screen_by_ligand_contact = False

for i in range(3):
    universe.run(snapshot_idx=i)





universe.screen_pockets()

pockets = []
for i in range(3):
    pockets.extend(universe.pockets(i))

for pocket in pockets:
    print(pocket.lining_atom_idx)







