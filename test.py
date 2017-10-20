import mdtraj
from AS import AS_Universe
from multiprocessing import Pool, cpu_count
import matplotlib
from timeit import default_timer

test_ligand_path = '/Users/haotian/Dropbox/pycharm_project/AlphaSpace/Test_system/lig.pdb'
test_protein_path = '/Users/haotian/Dropbox/pycharm_project/AlphaSpace/Test_system/prot.pdb'


traj_path = '/Users/haotian/Dropbox/Structure/Mdm2/prod_100ns_1000ss.mdcrd'
top_path = '/Users/haotian/Dropbox/Structure/Mdm2/1YCR_solv.prmtop'

#
# start = default_timer()
# traj = mdtraj.load(traj_path,top=top_path,stride = 100)
# print(start - default_timer())
#
# universe = AS_Universe(receptor=traj)




start = default_timer()
lig = mdtraj.load(test_ligand_path)
prot = mdtraj.load(test_protein_path)

universe = AS_Universe(receptor=prot,binder=lig)

universe.run()

print()
