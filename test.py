import mdtraj
from AS import *
from multiprocessing import Pool, cpu_count
import matplotlib
from timeit import default_timer

test_ligand_path = '/Users/haotian/Dropbox/pycharm_project/AlphaSpace/Test_system/lig.pdb'
test_protein_path = '/Users/haotian/Dropbox/pycharm_project/AlphaSpace/Test_system/prot.pdb'


traj_path = '/Users/haotian/Dropbox/Structure/Mdm2/prod_100ns_1000ss.mdcrd'
top_path = '/Users/haotian/Dropbox/Structure/Mdm2/1YCR_solv.prmtop'


start = default_timer()
traj = mdtraj.load(traj_path,top=top_path,stride = 100)
print(start - default_timer())

session = AS_Session(receptor=traj)


#
# start = default_timer()
# traj = mdtraj.load(traj_path,top=top_path,stride = 100)
# print(start - default_timer())
# session = AS_Session(receptor=traj)

# session = AS_Session(mdtraj.load(test_protein_path),mdtraj.load(test_ligand_path))

for i in range(session.n_frames):
    session.run()
