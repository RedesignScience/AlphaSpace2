import numpy as np

import scipy.cluster.hierarchy as hier
import scipy.sparse.csgraph as csg
import scipy.spatial.distance as scidist

import alphaspace2.Features as features

from .View import gen_pdb_line

import secrets
import string


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
#            print(p,np.sum(space_dict[p]))
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

def extendedCommunityAnalysis(protein, snapshot, atoms_to_writeout, save_npy, output_prefix):
  
 
# This function wraps community analysis into a single call. It returns a tuple of three elements. The first is a list of community properties in the form of
# nonploar space, volume, occluded ASA for community 1, nonploar space, volume, occluded ASA for community 2, etc. The second is a dictionary containing the pockets 
# making up the individual communities; the third is a list of all the pockets in the snapshot. The function can be used to write out all the alfa and/or beta atoms
# in the communities in pdb format with each community as a separate residue, containing the nonpolar space of the community in the b factor columns. This makes
# visualization easy and convenient. The user can also have the community properties saved as an npy file for subsequent retrieval and usage by specifying 
# save_npy = True in the function call. See the corresponding cookbook for usage examples.
    
    community_properties = []
    
    temp_space_dict = {}
    temp_coords_dict = {}
    temp_contact_dict = {}
    temp_score_dict = {}
    temp_centroids_dict = {}
    for px,pocket in enumerate(snapshot.pockets):
        temp_coords_dict[px] = [list(b.xyz) for b in pocket.betas]
        temp_score_dict[px] = np.array([0 for b in  pocket.betas]) 
        temp_space_dict[px] = np.array([b.nonpolar_space for b in  pocket.betas]) # ANALYSIS IS BASED ON NONPOLAR SPACE
        temp_contact_dict[px] = list(pocket.lining_atoms_idx)
        temp_centroids_dict[px] = pocket.centroid
    prot_coords = protein.xyz[0]*10   ### change nm to Angstrom

    surface_communities = genCommunityPocket(prot_coords, temp_centroids_dict, temp_space_dict, temp_contact_dict, temp_score_dict, corecut = 100, auxcut = 30, tight_option = True, tight_cutoff_core = 12.5, tight_cutoff_aux = 6.5)

    surface_communities_props = {}
    
    for cx,community in surface_communities.items():
        temp_coords_array = []
        for pock in community['core_pockets'] + community['aux_pockets']:
            temp_coords_array.extend(temp_coords_dict[pock])
            
        temp_coords_array = np.array(temp_coords_array)
        volume = features._get_grid_volume(temp_coords_array)
        occluded_asa = features._get_pharmacophore_fingerprint(protein,temp_coords_array)
        surface_communities_props[cx] = {}
        surface_communities_props[cx]['space'] = community['space']
        
        community_properties.append(round(community['space']))
        community_properties.append(round(volume))
        community_properties.append(round(occluded_asa['Total_OASA']))
        
    all_pockets = list(snapshot.pockets)
    
    if save_npy == True:
        np.save(output_prefix + '.npy', community_properties)
        
    if atoms_to_writeout == 'alfas':

        # handle is opened in append mode in case we want to place the output from multiple snapshots into a single, multimodel PDB file; hence the MODEL
        # and ENDMODEL lines.

        handle = open(output_prefix + '_community_alfas.pdb', 'a')
        handle.write("MODEL\n")

        atom_index = 1

        for i in surface_communities:
            
            core_pockets = list(surface_communities[i]['core_pockets'])
            aux_pockets = list(surface_communities[i]['aux_pockets'])

            # atom names are random, 4-character alfanumeric strings (lowercase letters included). This is done because having atoms with identical names and
            # different b factors can lead to erratic behavior with some visualization tools.
            for j in core_pockets:
                pocket = all_pockets[j]
                alfas = list(pocket.alphas)
                for alfa in alfas:
                    handle.write(gen_pdb_line(atomIndex=atom_index, atomName=''.join(secrets.choice(string.ascii_uppercase + string.ascii_lowercase + string.digits) for i in range(4)), resName='ACP', resIndex=i+1, chainName=" ", bfactor=surface_communities[i]['space'], element="C", xyz=alfa.xyz, occupancy=alfa.nonpolar_space))
                    atom_index += 1
                                            
            for j in aux_pockets:
                pocket = all_pockets[j]
                alfas = list(pocket.alphas)
                for alfa in alfas:
                    handle.write(gen_pdb_line(atomIndex=atom_index, atomName=''.join(secrets.choice(string.ascii_uppercase + string.ascii_lowercase + string.digits) for i in range(4)), resName='AAP', resIndex=i+1, chainName=" ", bfactor=surface_communities[i]['space'], element="C", xyz=alfa.xyz, occupancy=alfa.nonpolar_space))
                    atom_index += 1
                        
        handle.write("ENDMDL\n")
        handle.close()
        
        return community_properties, surface_communities, all_pockets

    elif atoms_to_writeout == 'betas':

        handle = open(output_prefix + '_community_betas.pdb', 'a')
        handle.write("MODEL\n")

        atom_index = 1
        
        for i in surface_communities:
            
            core_pockets = list(surface_communities[i]['core_pockets'])
            aux_pockets = list(surface_communities[i]['aux_pockets'])
                        
            for j in core_pockets:
                pocket = all_pockets[j]
                betas = list(pocket.betas)
                for beta in betas:
                    handle.write(gen_pdb_line(atomIndex=atom_index, atomName=''.join(secrets.choice(string.ascii_uppercase + string.ascii_lowercase + string.digits) for i in range(4)), resName='BCP', resIndex=i+1, chainName=" ", bfactor=surface_communities[i]['space'], element="C", xyz=beta.xyz, occupancy=beta.nonpolar_space))
                    atom_index += 1
                        
            for j in aux_pockets:
                pocket = all_pockets[j]
                betas = list(pocket.betas)
                for beta in betas:
                    handle.write(gen_pdb_line(atomIndex=atom_index, atomName=''.join(secrets.choice(string.ascii_uppercase + string.ascii_lowercase + string.digits) for i in range(4)), resName='BAP', resIndex=i+1, chainName=" ", bfactor=surface_communities[i]['space'], element="C", xyz=beta.xyz, occupancy=beta.nonpolar_space))
                    atom_index += 1

        handle.write("ENDMDL\n")
        handle.close()
        
        return community_properties, surface_communities, all_pockets
    
    elif atoms_to_writeout == 'both':
        
        handle = open(output_prefix + '_community_betas.pdb', 'a')
        handle.write("MODEL\n")

        atom_index = 1
        
        for i in surface_communities:
            
            core_pockets = list(surface_communities[i]['core_pockets'])
            aux_pockets = list(surface_communities[i]['aux_pockets'])

            for j in core_pockets:
                pocket = all_pockets[j]
                betas = list(pocket.betas)
                for beta in betas:
                    handle.write(gen_pdb_line(atomIndex=atom_index, atomName=''.join(secrets.choice(string.ascii_uppercase + string.ascii_lowercase + string.digits) for i in range(4)), resName='BCP', resIndex=i+1, chainName=" ", bfactor=surface_communities[i]['space'], element="C", xyz=beta.xyz, occupancy=beta.nonpolar_space))
                    atom_index += 1
                        
            for j in aux_pockets:
                pocket = all_pockets[j]
                betas = list(pocket.betas)
                for beta in betas:
                    handle.write(gen_pdb_line(atomIndex=atom_index, atomName=''.join(secrets.choice(string.ascii_uppercase + string.ascii_lowercase + string.digits) for i in range(4)), resName='BAP', resIndex=i+1, chainName=" ", bfactor=surface_communities[i]['space'], element="C", xyz=beta.xyz, occupancy=beta.nonpolar_space))
                    atom_index += 1
                
        handle.write("ENDMDL\n")
        handle.close()
    
        handle = open(output_prefix + '_community_alfas.pdb', 'a')
        handle.write("MODEL\n")
        
        atom_index = 1

        for i in surface_communities:

            core_pockets = list(surface_communities[i]['core_pockets'])
            aux_pockets = list(surface_communities[i]['aux_pockets'])
            
            for j in core_pockets:
                pocket = all_pockets[j]
                alfas = list(pocket.alphas)
                for alfa in alfas:
                    handle.write(gen_pdb_line(atomIndex=atom_index, atomName=''.join(secrets.choice(string.ascii_uppercase + string.ascii_lowercase + string.digits) for i in range(4)), resName='ACP', resIndex=i+1, chainName=" ", bfactor=surface_communities[i]['space'], element="C", xyz=alfa.xyz, occupancy=alfa.nonpolar_space))
                    atom_index += 1
                            
            for j in aux_pockets:
                pocket = all_pockets[j]
                alfas = list(pocket.alphas)
                for alfa in alfas:
                    handle.write(gen_pdb_line(atomIndex=atom_index, atomName=''.join(secrets.choice(string.ascii_uppercase + string.ascii_lowercase + string.digits) for i in range(4)), resName='AAP', resIndex=i+1, chainName=" ", bfactor=surface_communities[i]['space'], element="C", xyz=alfa.xyz, occupancy=alfa.nonpolar_space))
                    atom_index += 1
                
        handle.write("ENDMDL\n")
        handle.close()
        
        return community_properties, surface_communities, all_pockets
    
    elif atoms_to_writeout == 'None':
        return community_properties, surface_communities, all_pockets

    else:
        raise ValueError("Incorrect output atoms specified. Possible options are alfas, betas both, or None.")

