import alphaspace
import numpy as np
import mdtraj


testProteins_bcl2 = ['/Users/haotian/Dropbox/pycharm_project/AlphaSpace/examples/test/bcl2/prot/{}.pdb'.format(i) for i in range(1, 1 + 10)]

testLigand_bcl2 = ['/Users/haotian/Dropbox/pycharm_project/AlphaSpace/examples/test/bcl2/lig/{}.pdb'.format(i) for i in range(1, 1 + 10)]

u = alphaspace.AS_Universe()
for _ in range(1):
    for test_protein_path, test_ligand_path in zip(testProteins_bcl2,testLigand_bcl2):
        # print(test_protein_path,test_ligand_path)
        ligand = mdtraj.load(test_ligand_path)
        protein = mdtraj.load(test_protein_path)
        u.set_receptor(protein, append=True)
        u.set_binder(ligand, append=True)
u.run_alphaspace()


abs_sasa = []
for si in u.snapshots_indices:
    beta_coords = []
    for pocket in u.pockets(si):
        for beta in pocket.betas:
            beta_coords.append(beta.centroid)

    covered_sasa = alphaspace.getSASA(u.receptor.traj[si],np.array(beta_coords))
    sasa = alphaspace.getSASA(u.receptor.traj[si])

    abs_sasa.append(sasa - covered_sasa)

    print(np.count_nonzero(np.any((sasa - covered_sasa) > 0)),flush=True)




#
# exposed_index =  np.arange(covered_sasa.shape[0])[abs_sasa > 0]
#
# print(exposed_index.shape)
