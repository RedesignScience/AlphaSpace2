import numpy as np

import scipy.cluster.hierarchy as hier
import scipy.sparse.csgraph as csg
import scipy.spatial.distance as scidist


def genCommunityPocket(prot_coords, centroid_dict, space_dict, contact_index_dict, score_dict, corecut=100, auxcut=30,
                       tight_option=True, tight_cutoff_core=12.5, tight_cutoff_aux=6.5):
    """
    Generate Pocket Communities one pocket snapshot
    :param prot_coords: coordinates of proteins from mdtraj object
    :param dicts of properties
    :param corecut: cutoff score for core pocket
    :param auxcut: cutoff score for  auxiliary pocket
    """
    pkt_list = space_dict.keys()
    pkt_core_list = []
    ftdict = {}
    for p in pkt_list:
        arr_set = contact_index_dict[p]
        prot_centr = np.mean(prot_coords[arr_set], axis=0)
        ftdict[p] = [arr_set, prot_centr, centroid_dict[p], [p], np.sum(space_dict[p]), np.sum(score_dict[p])]
        if np.sum(space_dict[p]) >= corecut:
            pkt_core_list.append(p)
    compkt = {}
    if len(pkt_core_list) > 1:
        coremat = []
        for ix, i in enumerate(pkt_core_list[:-1]):
            for jx, j in enumerate(pkt_core_list[(ix + 1):]):
                if np.any(np.intersect1d(ftdict[i][0], ftdict[j][0])) and \
                        (sum((ftdict[i][2] - ftdict[i][1]) * (ftdict[j][2] - ftdict[j][1])) > 0 or (
                                np.linalg.norm(ftdict[i][1] - ftdict[j][1]) > np.linalg.norm(
                            ftdict[i][2] - ftdict[j][2]))):
                    coremat.append(1)
                else:
                    coremat.append(0)
        coremat = scidist.squareform(np.array(coremat))
        nclust, clust = csg.connected_components(coremat, directed=False)
        if tight_option:
            # Do average linkage clustering
            nclust, clust = CoreCluster(pkt_core_list, ftdict, nclust, clust, CC_cut=tight_cutoff_core)
        for i in range(nclust):
            compkt[i] = [[], [], [], 0, 0]
        for i, j in zip(clust, pkt_core_list):
            compkt[i][0].append(j)
    elif len(pkt_core_list) == 1:
        nclust = 1
        compkt[0] = [[pkt_core_list[0]], [], [], 0, 0]
    elif len(pkt_core_list) == 0:
        nclust = 0
        print("No Pocket with score > %d!" % (corecut))
    for i in range(nclust):
        for j in compkt[i][0]:
            for k in pkt_list:
                if k not in compkt[i][0]:
                    if tight_option:
                        if np.any(np.intersect1d(ftdict[j][0], ftdict[k][0])) > 0 and \
                                (sum((ftdict[j][2] - ftdict[j][1]) * (ftdict[k][2] - ftdict[k][1])) > 0 or \
                                 (np.linalg.norm(ftdict[j][1] - ftdict[k][1]) > np.linalg.norm(
                                     ftdict[j][2] - ftdict[k][2]))) and \
                                np.linalg.norm(ftdict[j][2] - ftdict[k][2]) < tight_cutoff_aux:
                            if ftdict[k][4] >= auxcut and ftdict[k][4] < corecut:
                                if k not in compkt[i][1]:
                                    compkt[i][1].append(k)
                            elif k not in compkt[i][2] and ftdict[k][4] < corecut:
                                compkt[i][2].append(k)
                    else:
                        if np.any(np.intersect1d(ftdict[j][0], ftdict[k][0])) > 0 and \
                                (sum((ftdict[j][2] - ftdict[j][1]) * (ftdict[k][2] - ftdict[k][1])) > 0 or \
                                 (np.linalg.norm(ftdict[j][1] - ftdict[k][1]) > np.linalg.norm(
                                     ftdict[j][2] - ftdict[k][2]))):
                            if ftdict[k][4] >= auxcut and ftdict[k][4] < corecut:
                                if k not in compkt[i][1]:
                                    compkt[i][1].append(k)
                            elif k not in compkt[i][2] and ftdict[k][4] < corecut:
                                compkt[i][2].append(k)
    for i in range(nclust):
        for j in compkt[i][0]:
            compkt[i][3] += ftdict[j][4]
            compkt[i][4] += ftdict[j][5]
        for j in compkt[i][1]:
            compkt[i][3] += ftdict[j][4]
            compkt[i][4] += ftdict[j][5]
    sortedcompkt = sorted(compkt.items(), key=lambda x: x[1][3], reverse=True)
    pocket_community_dict = {}
    for ax, (a, b) in enumerate(sortedcompkt):
        temp_dict = {'core_pockets': b[0]}
        temp_dict['aux_pockets'] = b[1]
        temp_dict['minor_pockets'] = b[2]
        temp_dict['space'] = round(b[3])
        temp_dict['score'] = round(b[4], 2)
        pocket_community_dict[ax] = temp_dict
    return pocket_community_dict


def CoreCluster(pkt_core_list, ftdict, nclust, clust, CC_cut=8.5):
    """Core Cluster by average linkage using distance cutoff
    :param pkt_core_list: the core pocket id list
    :param ftdict: dict of pocket information
    :param nclust: number of original cluster
    :param clust: cluste index
    :param CC_cut: the cutoff distance
    """
    pktdict = dict(zip(pkt_core_list, clust))
    for i in range(nclust):
        tmp_core_list = [key for key, value in pktdict.items() if value == i]
        if len(tmp_core_list) > 1:
            tmp_core_center = [ftdict[j][2] for j in tmp_core_list]
            tmp_zmat = hier.linkage(tmp_core_center, method="average")
            tmp_cluster = hier.fcluster(tmp_zmat, CC_cut, criterion='distance')
            tmp_ncluster = max(tmp_cluster)
            if tmp_ncluster > 1:
                for m, n in zip(tmp_core_list, tmp_cluster):
                    if n > 1:
                        pktdict[m] = nclust + n - 2
                nclust = nclust + tmp_ncluster - 1
    newcluster = [pktdict[i] for i in pkt_core_list]

    return nclust, newcluster
