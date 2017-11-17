import mdtraj
import sys
from alphaspace import *


# for _ in range(10):
#     for i in range(1, 1 + 10):
#         test_ligand_path = './bcl2/lig/{}.pdb'.format(i)
#         test_protein_path = './bcl2/prot/{}.pdb'.format(i)
#         ligand = mdtraj.load(test_ligand_path)
#         protein = mdtraj.load(test_protein_path)
#         universe.set_receptor(protein, append=True)
#         universe.set_binder(ligand, append=True)

traj = mdtraj.load('/Users/haotian/Desktop/Nemo-MD/28-6/Prod_1/prod_1x100.mdcrd',
                   top='/Users/haotian/Desktop/Nemo-MD/28-6/3CL3_solv.prmtop')
universe = AS_Universe(traj)
print('loading finished')

universe.run_alphaspace_mp(4)
