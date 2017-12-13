import numpy as np
from scipy import spatial

# Import hp-type and load dictionary
with open('hp_types_dict.dat', 'r') as f:
    templines = f.readlines()

types_dict = {}
for t in templines[1:]:
    temp = t.split(',')
    if temp[0] in types_dict.keys():
        types_dict[temp[0]][temp[1]] = [temp[2].strip(), temp[3].strip()]
    else:
        types_dict[temp[0]] = {}
        types_dict[temp[0]][temp[1]] = [temp[2].strip(), temp[3].strip()]

# polar_atoms=np.array(['OA','N','O','S','NS'])
polar_atoms = np.array(['OA', 'OS', 'N', 'NS', 'NA', 'S', 'SA'])

# Import autodock_types_dict
autodock_types_dict = {}
with open('autodock_atom_type_info.dat', 'r') as f:
    templines = f.readlines()
for t in templines[1:]:
    temp = t.split(',')
    autodock_types_dict[temp[0].strip()] = [float(temp[1]) / 2.0, True if int(temp[7]) else False]

# Import vina parameters
with open('vina_params.dat', 'r') as f:
    templines = f.readlines()
vina_weights_dict = {}
for t in templines[1:]:
    temp = t.split(',')
    vina_weights_dict[temp[2].split()[0]] = float(temp[1])


def NP_interp(r):
    """
    #step function for nonpolar interactions.
    :param r: radii
    :type r: fload
    :return: stepped
    :rtype: float
    """
    if r < 0.5:
        x = 1.0
    elif r > 1.5:
        x = 0.0
    else:
        x = 1.5 - r
    #        x=np.interp(r, [0.i5,1.5], [1,0])
    return x


def P_interp(r):  ##step for polar
    if r < -0.7:
        x = 1.0
    elif r >= 0:
        x = 0.0
    else:
        x = -r / 0.7
    #        x=np.interp(r, [0.i5,1.5], [1,0])
    return x


def calc_score(g1, g2, rep, hydrophobe, hbond):
    """
    Sum different components of vina score

    :param g1: gaussian 1
    :param g2: gaussian 2
    :param rep: repulsion
    :param hydrophobe: hydrophobic
    :param hbond: hydrogen bond
    :return: summed score
    :rtype: float
    """
    return g1 * vina_weights_dict['gauss_1'] \
           + g2 * vina_weights_dict['gauss_2'] \
           + rep * vina_weights_dict['repulsion'] \
           + hydrophobe * vina_weights_dict['hydrophobic'] \
           + hbond * vina_weights_dict['Hydrogen']


def assign_hp(prot_coord, prot_types, hp_type):
    """
    This assigns hydrophobicity type for UNK atoms

    :param prot_coord: protein coordinate shape =  N*3
    :param prot_types: protein atom types shape = N
    :param hp_type: type_to be modified.
    :return: hp_type
    """

    tree = spatial.KDTree(prot_coord)
    for ix in np.where(hp_type == 'UNK')[0]:
        indx = np.array(tree.query_ball_point(prot_coord[ix], 2.0))
        if prot_types[ix] not in polar_atoms:
            if np.any(np.in1d(prot_types[indx[indx != ix]], polar_atoms)):
                hp_type[ix] = 'NNP'
            else:
                hp_type[ix] = 'NP'
        else:
            hp_type[ix] = 'XXX'
    return hp_type


def assign_acc(prot_coord, prot_types, hp_type):
    """

    :param prot_coord:
    :param prot_types:
    :param hp_type:
    :return:
    """
    tree = spatial.KDTree(prot_coord)
    for ix in np.where(hp_type == 'UNK')[0]:
        indx = np.array(tree.query_ball_point(prot_coord[ix], 2.0))
        if prot_types[ix] in ['OA', 'OS', 'SA', 'S']:
            hp_type[ix] = 'P'
        else:
            hp_type[ix] = 'XXX'
    return hp_type


def assign_don(prot_coord, prot_types, hp_type):
    tree = spatial.KDTree(prot_coord)
    for ix in np.where(hp_type == 'UNK')[0]:
        indx = np.array(tree.query_ball_point(prot_coord[ix], 2.0))
        if prot_types[ix] in ['N', 'NS', 'NA']:
            if len(indx[indx != ix]) > 2:
                hp_type[ix] = 'NPP'
            else:
                hp_type[ix] = 'P'
        else:
            hp_type[ix] = 'XXX'
    return hp_type


ali_atoms = np.array(['C', 'A'])


def pre_process_pdbqt(pdbqt_file):
    ali_atoms = np.array(['C', 'A'])
    with open(pdbqt_file, 'r') as f:
        templines = f.readlines()
    hp_type = []
    don_type = []
    acc_type = []
    prot_coord = []
    prot_types = []
    for t in templines:
        if t[0:6] in ['ATOM  ', 'HETATM'] and t[-3:].strip() not in ['H', 'HD', 'HS']:
            prot_coord.append([float(t[30:38].strip()), float(t[38:46].strip()), float(t[46:54].strip())])
            prot_types.append(t[76:].strip())
            if t[17:20].strip() in types_dict.keys():
                if t[12:16].strip() in types_dict[t[17:20].strip()].keys():
                    if t[76:].strip() in ali_atoms:
                        hp_type.append(types_dict[t[17:20].strip()][t[12:16].strip()][0])
                        don_type.append('XXX')
                        acc_type.append('XXX')
                    if t[76:].strip() in polar_atoms:
                        don_type.append(types_dict[t[17:20].strip()][t[12:16].strip()][0])
                        acc_type.append(types_dict[t[17:20].strip()][t[12:16].strip()][1])
                        hp_type.append('XXX')
                else:
                    hp_type.append('UNK')
                    don_type.append('UNK')
                    acc_type.append('UNK')
            else:
                hp_type.append('UNK')
                don_type.append('UNK')
                acc_type.append('UNK')
    prot_coord = np.array(prot_coord)
    prot_types = np.array(prot_types)
    hp_type = np.array(hp_type)
    acc_type = np.array(acc_type)
    don_type = np.array(don_type)
    if np.any(hp_type == 'UNK'):
        assign_hp(prot_coord, prot_types, hp_type)
        assign_don(prot_coord, prot_types, don_type)
        assign_acc(prot_coord, prot_types, acc_type)
    return prot_coord, prot_types, hp_type, acc_type, don_type


def prep_temp_dict(type_list):
    temp_prb_dict = {}
    for t in type_list:
        temp_prb_dict[t] = []
    return temp_prb_dict


nrdld_pscore_dict = {}

for i in ['train', 'val']:
    nrdld_pscore_dict[i] = {}
    for j in ['drug', 'nondrug']:
        nrdld_pscore_dict[i][j] = {}
        for k in nrdld_list_dict[i][j]:
            try:
                nrdld_pscore_dict[i][j][k] = {}
                prot_coord, prot_types, hp_type, acc_type, don_type = pre_process_pdbqt(
                    i + '/' + j + '/' + k + '/protein.pdbqt')  ### this part pdbqt
                for pxx in space_communities_dict[i][j][k].keys():
                    pock_list = [pp for pp in
                                 space_communities_dict[i][j][k][pxx][0] + space_communities_dict[i][j][k][pxx][1] +
                                 space_communities_dict[i][j][k][pxx][2]]  ## where to feed coordinates
                    coords_list = []
                    for pp in pock_list:
                        coords_list.extend(beta_space_coords_dict[i][j][k][pp][0])
                    dist = scipy.spatial.distance.cdist(coords_list, prot_coord)
                    temp_prb_dict = prep_temp_dict(['C', 'Br', 'F', 'Cl', 'I', 'OA', 'SA', 'N', 'P'])
                    for px in range(0, len(coords_list)):
                        dist_bool = dist[px] <= 8.0
                        temp_dist = dist[px][dist_bool]
                        temp_type = prot_types[dist_bool]
                        NP_type = (hp_type[dist_bool] == 'NP')
                        Pdon_type = (don_type[dist_bool] == 'P')
                        Pacc_type = (acc_type[dist_bool] == 'P')
                        dist_radii = np.array([autodock_types_dict[ty][0] for ty in temp_type])
                        for prb in ['C', 'Br', 'F', 'Cl', 'I', 'OA', 'SA', 'N', 'P']:
                            probe_dist = autodock_types_dict[prb][0]
                            proc_dist = temp_dist - dist_radii - probe_dist
                            g1 = np.sum(np.exp(-(proc_dist / 0.5) ** 2))
                            g2 = np.sum(np.exp(-((proc_dist - 3.0) / 2.0) ** 2))
                            rep = np.sum([dd ** 2 if dd < 0.0 else 0.0 for dd in proc_dist])
                            if prb in ['C', 'Br', 'Cl', 'F', 'I']:
                                h1 = np.sum([NP_interp(dd) for dd in proc_dist[NP_type]])
                                h2 = 0.0
                            elif prb in ['OA', 'SA']:
                                h1 = 0.0
                                h2 = np.sum([P_interp(dd) for dd in proc_dist[Pdon_type]])
                            elif prb in ['N', 'P']:
                                h1 = 0.0
                                h2 = np.sum([P_interp(dd) for dd in proc_dist[Pacc_type]])
                            temp_prb_dict[prb].append([calc_score(g1, g2, rep, h1, h2), g1, g2, rep, h1, h2])
                    nrdld_pscore_dict[i][j][k][pxx] = temp_prb_dict
            except:
                print('ERR', i, j, k)
            print(i, j, k)
