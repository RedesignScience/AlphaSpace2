



def get_asa(trajectory,):

    beta_coords = []
    for pocket in u.pockets(si):
        for beta in pocket.betas:
            beta_coords.append(beta.centroid)

    covered_sasa = alphaspace.getSASA(u.receptor.traj[si], np.array(beta_coords))
    sasa = alphaspace.getSASA(u.receptor.traj[si])

    abs_sasa.append(sasa - covered_sasa)

    print(np.count_nonzero(np.any((sasa - covered_sasa) > 0)), flush=True)





