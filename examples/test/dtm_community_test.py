import alphaspace as al
import matplotlib.pyplot as plt
import numpy as np

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
    u._gen_community(100)

    # Plotting similarities.
    # sim1 = []
    # sim2 = []
    #
    # for snapshot_idx in u.snapshots_indices:
    #     c_sims = [c.similarity_atoms(lining_atoms) for c in u._communities[snapshot_idx]]
    #     c_sims.sort(reverse=True)
    #     if len(c_sims) > 1:
    #         sim1.append(c_sims[0])
    #         sim2.append(c_sims[1])
    #     else:
    #         sim1.append(c_sims[0])
    #         sim2.append(0)
    #
    # plt.ylim([0,0.5])
    # plt.plot(range(160), sim1,label = '1st community')
    # plt.plot(range(160), sim2,label = '2nd community')
    # legend = plt.legend(loc='upper right', shadow=True)
    # plt.title("{} {} community similarity to binding community".format(protein, system))
    # plt.savefig('/Users/haotian/Desktop/plots/{}.{}.png'.format(protein,system),tight = True)
    # print(protein,system)

    # Plotting aggregated pocket score/space

    # sim1_score = []
    # sim2_score = []
    #
    # for snapshot_idx in u.snapshots_indices:
    #     sim_sort_com = sorted([c for c in u._communities[snapshot_idx]],key = lambda x: x.similarity_atoms(lining_atoms),reverse=True)
    #
    #     if len(sim_sort_com) > 1:
    #         sim1_score.append(np.sum([p.space for p in sim_sort_com[0].pockets]))
    #         sim2_score.append(np.sum([p.space for p in sim_sort_com[1].pockets]))
    #     else:
    #         sim1_score.append(np.sum([p.space for p in sim_sort_com[0].pockets]))
    #         sim2_score.append(0)
    #
    #
    # # plt.ylim([0, 0.5])
    # plt.plot(range(160), sim1_score, label='1st community')
    # plt.plot(range(160), sim2_score, label='2nd community')
    # legend = plt.legend(loc='upper right', shadow=True)
    # plt.title("{} {} agg space for 1st 2nd similiar communities".format(protein, system))
    # plt.savefig('/Users/haotian/Desktop/plots//space/{}.{}.png'.format(protein, system), tight=True)
    # print(protein, system)

    # Plot aggregated score of best scored and most similar community

    # sim1_sim = []
    # score1_sim = []
    #
    # for snapshot_idx in u.snapshots_indices:
    #     sim_sort_com = sorted([c for c in u._communities[snapshot_idx]],key = lambda x: x.similarity_atoms(lining_atoms),reverse=True)
    #     score_sort_com = sorted([c for c in u._communities[snapshot_idx]],key = lambda x: np.sum([p.score for p in x.pockets]),reverse=False)
    #
    #
    #     # if len(sim_sort_com) > 1:
    #     sim1_sim.append(sim_sort_com[0].similarity_atoms(lining_atoms))
    #     score1_sim.append(score_sort_com[0].similarity_atoms(lining_atoms))
    #     # else:
    #     #     sim1_score.append(np.sum([p.space for p in sim_sort_com[0].pockets]))
    #     #     sim2_score.append(0)
    #
    #
    # # plt.ylim([0, 0.5])
    # plt.plot(range(160), sim1_sim, label='Most similar community')
    # plt.plot(range(160), score1_sim, label='Best Scoring community')
    # legend = plt.legend(loc='upper right', shadow=True)
    # plt.ylabel('Tanimoto similarity to crystal binding community')
    # plt.title("{} {} Best Scoring community to Most Similar community".format(protein, system))
    # plt.savefig('/Users/haotian/Desktop/plots/scored_sim/{}.{}.png'.format(protein, system), tight=True)
    # print(protein, system)

    # Score vs similarity

    # communities = []
    #
    # for snapshot_idx in u.snapshots_indices:
    #     communities += [c for c in u._communities[snapshot_idx]]
    #
    # plt.scatter([c.similarity_atoms(lining_atoms) for c in communities],[np.sum([p.score for p in c.pockets]) for c in communities])
    #
    # # legend = plt.legend(loc='upper right', shadow=True)
    # plt.xlabel('Tanimoto similarity to crystal binding community')
    # plt.ylabel('Aggregated Score')
    # plt.title("{} {} Score vs Similarity".format(protein, system))
    # plt.savefig('/Users/haotian/Desktop/plots/score_vs_sim/{}.{}.png'.format(protein, system), tight=True)
    # print(protein, system)


    # all attributes of community in interface
    communities = []

    for snapshot_idx in u.snapshots_indices:
        communities += [c for c in u._communities[snapshot_idx]]

    sim1_score = []
    sim1_space = []

    for snapshot_idx in u.snapshots_indices:
        sim_sort_com = max([c for c in u._communities[snapshot_idx]],key = lambda x: x.similarity_atoms(lining_atoms))
        # score_sort_com = sorted([c for c in u._communities[snapshot_idx]],key = lambda x: np.sum([p.score for p in x.pockets]),reverse=False)

        sim1_score.append(np.sum([p.score for p in sim_sort_com.pockets]))
        sim1_space.append(np.sum([p.space for p in sim_sort_com.pockets]))

    score_average = np.average(sim1_score)
    score_std = np.std(sim1_score)
    space_average = np.average(sim1_space)
    space_std = np.std(sim1_space)


    print(protein,system,score_average,score_std,space_average,space_std)