import mdtraj
from AS import *
from multiprocessing import Pool, cpu_count
import matplotlib
from timeit import default_timer

test_ligand_path = '/Users/haotian/Dropbox/pycharm_project/AlphaSpace/Test_system/lig.pdb'
test_protein_path = '/Users/haotian/Dropbox/pycharm_project/AlphaSpace/Test_system/prot.pdb'

lig = mdtraj.load(test_ligand_path)
prot = mdtraj.load(test_protein_path)

session = AS_Session(prot, lig)

session.run()
# session.screen_by_ligand_contact()




# session._run(cpu_count())
#
# cluster = session.receptor.clusters[0]
#
# cluster.screen_by_contact()
