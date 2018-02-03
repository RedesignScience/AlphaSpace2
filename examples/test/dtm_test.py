import alphaspace as al
import matplotlib.pyplot as plt

index = (('bcl2', 'apo1'), ('bcl2', 'apo2'), ('bcl2', 'holo1'), ('bclxl', 'apo1'), ('bclxl', 'apo2'), ('bclxl', 'apo3'),
         ('bclxl', 'apo4'), ('bclxl', 'apo5'), ('bclxl', 'holo1'), ('bclxl', 'holo2'), ('bclxl', 'holo3'),
         ('bclxl', 'holo4'), ('hiv', 'apo1'), ('hiv', 'apo2'), ('hiv', 'holo1'), ('hpv', 'apo1'), ('hpv', 'holo1'),
         ('il2', 'apo1'), ('il2', 'apo2'), ('il2', 'holo1'), ('mdm4', 'apo1'), ('mdm4', 'holo1'), ('menin', 'apo1'),
         ('menin', 'apo2'), ('menin', 'holo1'), ('tnfa', 'apo1'), ('tnfa', 'holo1'), ('tnfr1', 'apo1'),
         ('tnfr1', 'apo2'), ('tnfr1', 'holo1'), ('xiap', 'apo1'), ('xiap', 'apo2'), ('xiap', 'apo3'), ('xiap', 'holo1'),
         ('xiap', 'holo2'), ('zipa', 'apo1'), ('zipa', 'apo2'), ('zipa', 'holo1'))

for protein, system in index:
    plt.clf()
    lining_atom_path = '/Users/haotian/Desktop/crystal_structures/{}.lining_atoms.dat'.format(protein)
    with open(lining_atom_path, 'r') as handle:
        lines = handle.readlines()
    lining_atoms = [int(i) for i in lines[0].split()]

    u = al.AS_Universe()
    u.load("/Users/haotian/Desktop/snapshots/{}/{}/vina_scored.as".format(protein, system))
    al.AS_Funct._gen_community(u)

    sim1 = []
    sim2 = []

    for snapshot_idx in u.snapshots_indices:
        c_sims = [c.similarity_atoms(lining_atoms) for c in u._communities[snapshot_idx]]
        c_sims.sort(reverse=True)
        if len(c_sims) > 1:
            sim1.append(c_sims[0])
            sim2.append(c_sims[1])
        else:
            sim1.append(c_sims[0])
            sim2.append(0)

    plt.ylim([0,0.5])
    plt.plot(range(160), sim1,label = '1st community')
    plt.plot(range(160), sim2,label = '2nd community')
    legend = plt.legend(loc='upper right', shadow=True)
    plt.title("{} {} community similarity to binding community".format(protein, system))
    plt.savefig('/Users/haotian/Desktop/plots/{}.{}.png'.format(protein,system),tight = True)
    print(protein,system)
