import alphaspace
import mdtraj
import _pickle as pickle




#
# universe = alphaspace.AS_Universe()
# for i in range(1, 1 + 10):
#     test_ligand_path = './tests/bcl2/lig/{}.pdb'.format(i)
#     test_protein_path = './tests/bcl2/prot/{}.pdb'.format(i)
#     ligand = mdtraj.load(test_ligand_path)
#     protein = mdtraj.load(test_protein_path)
#
#     universe.set_receptor(protein, append=True)
#     universe.set_binder(ligand, append=True)
# universe.run_alphaspace()

with open("/Users/haotian/Desktop/il2/apo1/md1/out.as",'rb') as handle:
    universe = pickle.load(handle)

for dp in universe.d_pockets:
    print(dp)


