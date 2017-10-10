import mdtraj
from AS import *
from multiprocessing import Pool, cpu_count
import matplotlib
from timeit import default_timer



traj_path = '/Users/haotian/Dropbox/Structure/Mdm2/prod_100ns_1000ss.mdcrd'
top_path = '/Users/haotian/Dropbox/Structure/Mdm2/1YCR_solv.prmtop'


start = default_timer()
traj = mdtraj.load(traj_path,top=top_path,stride = 100)
print(start - default_timer())

session = AS_Session(receptor=traj)

for i in range(session.n_frames):
    session.run()


# make a list of all pockets
pocket_list = []
for ss_idx in range(session.n_frames):
    for pocket in session.get_pockets(ss_idx):
        pocket_list.append(pocket)





# calculate pocket difference based on surface atom composition

# cluster the pockets based on surface atom overlap