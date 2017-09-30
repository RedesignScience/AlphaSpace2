import mdtraj
from AS import *
from multiprocessing import Pool,cpu_count


test_ligand_path = '/Users/haotian/Dropbox/pycharm_project/AlphaSpace/TestPDB/lig.pdb'
test_protein_path = '/Users/haotian/Dropbox/pycharm_project/AlphaSpace/TestPDB/prot.pdb'

lig = mdtraj.load(test_ligand_path)
prot = mdtraj.load(test_protein_path)


session = AS_Session(prot,lig)

session.run()

session.screen_by_ligand_contact()

for pocket in session.get_pockets():
    print(pocket._total_score*1000)

# session._run(cpu_count())
#
# cluster = session.receptor.clusters[0]
#
# cluster.screen_by_contact()




