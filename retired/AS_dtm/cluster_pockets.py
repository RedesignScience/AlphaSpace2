import pickle
from itertools import combinations

import matplotlib.pyplot as plt
import mdtraj as md
import networkx as nx
from scipy.spatial.distance import *

from DTM_class import *
from retired.AlphaSpace import *


def parse_traj(top_file,traj_file,num_snapshot):
    trajectory = md.load(traj_file,top=top_file)
    trajectory.remove_solvent(inplace=True)
    print(trajectory.n_chains)
    sliced_traj = trajectory[::int(trajectory.n_frames / num_snapshot)]
    os.chdir('tmp_as')
    pdb_list = []
    for i,snapshot in enumerate(sliced_traj):
        snapshot.save_pdb('{}.pdb'.format(i),)
        pdb_list.append('{}.pdb'.format(i))
    pdb_to_as_pickle(pdb_list=pdb_list,pickle_file='../as.pickle')
    os.chdir('..')
    # shutil.rmtree('tmp_as')
    return


def get_cluster(dist_matrix,p,method="average"):
    def fill_p_verts(p_verts,cluster):
        """
        Function: fill an empty list of lists (length equal to the number of clusters) with vertices_coords indices according to
        the clustering (becomes centroids in super clustering) Parameters p_verts: empty list of lists cluster: cluster
        indices
        """

        for i,clust_idx in enumerate(cluster):
            p_idx = clust_idx - 1
            p_verts[p_idx].append(i)
        return

    zmat = hier.linkage(dist_matrix,method=method)
    cluster = hier.fcluster(zmat,p,criterion='distance')
    clusters = [[] for _ in range(max(cluster))]
    fill_p_verts(clusters,cluster)
    return clusters


def pdb_to_as_pickle(pdb_list,pickle_file='as.pickle',param_file="AS_param.txt",option_file="AS_option.txt",):
    """
    This function reads in pdb file and process it with AlphaSpace, then save the needed data into a pickle object, you
    can retrieve it with the load_pickle function
    :param pdb_list: list of pdb files to be processed
    :param pickle_file: output pickle file location
    :param param_file: for running AlphaSpace
    :param option_file: for running AlphaSpace
    :return: None
    """

    output_object = []
    for pdb_file in pdb_list:
        with open(pdb_file,'r') as handle:
            pdblinelist = handle.readlines()
        os.system('mkdir {}'.format(pdb_file + '_as'))
        os.chdir(format(pdb_file + '_as'))
        pdb,ss = runAlphaSpace(pdblinelist,param_file=param_file,option_file=option_file,default=False,
                               isWrite=False)
        output_object.append([pdb,ss])
        os.chdir('..')
    with open(pickle_file,'wb') as handle:
        pickle.dump(output_object,handle)

    return


def load_as_pickle(pickle_file='as.pickle'):
    """
    This function simply load the AlphaSpace pickle object, and return a trajectory instance which contains all the
    information needed.
    :param pickle_file: pickle file location
    :return: AlphaSpace object list
    """
    with open(pickle_file,'rb') as handle:
        as_list = pickle.load(handle)
    return as_list


def cluster_pocket(as_list):
    master_pocket_index = []  # snapshot idx, pocket idx
    master_dpocket_list = []
    master_pocket_atom_list = []
    for snapshot_idx,snapshot_as in enumerate(as_list):
        pdb,snapshot = snapshot_as
        pockets = snapshot.sub_pockets
        for pocket in pockets:
            pocket.ss_idx = snapshot_idx
            pocket.ss = snapshot
            pocket_atom_list = np.zeros(len(pdb.prot))
            for atom in pocket.atoms:
                pocket_atom_list[atom] = 1
            master_pocket_atom_list.append(pocket_atom_list)
        master_pocket_index.extend([(snapshot_idx,_) for _ in range(len(pockets))])
    master_pocket_atom_list = np.array(master_pocket_atom_list)
    pocket_dist = pdist(master_pocket_atom_list,'jaccard')
    zmat = hier.linkage(pocket_dist,method="average")
    clust_dist = 0.75 * max(zmat[:,2])
    cluster = hier.fcluster(zmat,clust_dist,criterion='distance')
    clusters = [[] for _ in range(max(cluster))]
    fill_p_verts(clusters,cluster)
    for cluster in clusters:
        sPocket_list = []
        for snapshot_idx,pocket_idx in [master_pocket_index[i] for i in cluster]:
            sPocket_list.append(as_list[snapshot_idx][1].sub_pockets[pocket_idx])
        master_dpocket_list.append(dPocket(sPocket_list))
    # refinement by proximate
    dpocket_dist_matrix = np.zeros((len(master_dpocket_list),len(master_dpocket_list)))
    for i,j in combinations(range(len(master_dpocket_list)),2):
        dpocket_dist_matrix[i][j] = dpocket_dist_matrix[j][i] = master_dpocket_list[i] - master_dpocket_list[j]
    clustered_dpocket_list = [[master_dpocket_list[i] for i in cluster_idx] for cluster_idx in
                              get_cluster(squareform(dpocket_dist_matrix),p=5)]
    return clustered_dpocket_list


def plot_correlation(clustered_pockets):
    refined_dpocket_list = []
    for dpockets in clustered_pockets:
        dp = dpockets[0]
        for i in dpockets[1:]:
            dp.merge(i)
        if len(dp) > 1:
            refined_dpocket_list.append(dp)
    refined_dpocket_list.sort(key=lambda _: len(_),reverse=True)
    print([len(i) for i in refined_dpocket_list])
    correlation_matrix = np.zeros((len(refined_dpocket_list),len(refined_dpocket_list)))
    for i in range(len(refined_dpocket_list)):
        for j in range(i,len(refined_dpocket_list)):
            correlation_matrix[i,j] = correlation_matrix[j,i] = len(
                    refined_dpocket_list[i].get_ss_idx().intersection(refined_dpocket_list[j].get_ss_idx())) ** 2 / (
                                                                    len(refined_dpocket_list[1].get_ss_idx()) * len(
                                                                            refined_dpocket_list[j].get_ss_idx()))

    fig,ax = plt.subplots()
    heatmap = ax.pcolor(correlation_matrix,cmap=plt.cm.Blues,alpha=0.8)
    plt.colorbar(heatmap)
    plt.show()


def get_patches(dpocket_list,core_score_cutoff=100):
    dpocket_network = nx.Graph()
    core_list = []
    minor_list = []
    for i,dp in enumerate(dpocket_list):
        dpocket_network.add_node(i,population=len(dp),scores=dp.get_scores(),pocket_type=None)
    dpocket_atom_list = []
    for dpocket in dpocket_list:
        surface_atom_list = set()
        for pocket in dpocket.get_sPockets():
            pocket_atom = set(pocket.atoms)
            surface_atom_list = surface_atom_list | pocket_atom
        dpocket_atom_list.append(surface_atom_list)
    for i,j in combinations(range(len(dpocket_list)),2):
        if len(dpocket_atom_list[i] & dpocket_atom_list[j]) > 0:
            dpocket_network.add_edge(i,j)
    for node in dpocket_network.nodes_iter(data=True):
        # noinspection PyTypeChecker
        if len([score for score in node[1]['scores'] if score >= core_score_cutoff]) > len(node[1]["scores"])*0.3:
            node[1]['pocket_type'] = 'core'
            core_list.append(node[0])
        else:
            node[1]['pocket_type'] = 'minor'
            minor_list.append(node[0])
    patch_dict = {key: [] for key in core_list}
    for minor_node in minor_list:
        closest_cores = [core_list[0]]
        closest_core_dist = nx.algorithms.shortest_path_length(dpocket_network,minor_node,core_list[0])
        for core in core_list[1:]:
            dist = nx.algorithms.shortest_path_length(dpocket_network,minor_node,core)
            if dist < closest_core_dist:
                closest_core_dist = dist
                closest_cores = [core]
            elif dist == closest_core_dist:
                closest_cores.append(core)
            else:
                pass
        for core in closest_cores:
            patch_dict[core].append(minor_node)
    return patch_dict


def check_contact(dpocket):
    from scipy.spatial.distance import cdist

    for pocket in dpocket.get_sPockets():
        coord_lig = np.array([pdbinfo.coord(l) for l in pocket.ss.lig_pdb])
        pocket_center = np.array(pocket.vert_cntr())
        if min(cdist([pocket_center],coord_lig)[0]) < 5.0:
            return True
    return False


def get_contact_atoms(pocket_list):
    contact_pocket_atoms = set()
    for dpocket in pocket_list:
        if check_contact(dpocket):
            for pocket in dpocket.get_sPockets():
                contact_pocket_atoms = contact_pocket_atoms | set(pocket.atoms)
            #             print(len(contact_pocket_atoms))
    return contact_pocket_atoms


if __name__ == '__main__':
    # parse_traj(sys.argv[1], sys.argv[2],10)
    # pdb_to_as_pickle(sys.argv[1:])
    as_list = load_as_pickle(sys.argv[1])
    print("load complete")
    clustered_pockets = cluster_pocket(as_list)
    print("cluster complete")
    plot_correlation(clustered_pockets=clustered_pockets)
