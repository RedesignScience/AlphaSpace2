import timeit
from itertools import combinations

import numpy as np
import scipy.cluster.hierarchy as hier
from scipy.spatial.distance import cdist,pdist,squareform

from ..PDB import *
from ..PDB.alphaspace import AlphaSpace

"""
Import this file if clustering is needed.

python structure_clustering folders
to cluster given folder structures
"""


def contact_detection(atom_list1, atom_list2, threshold=1.3):
    """
    use numpy where and scipy.spatial.distance.cdist for fast process
    """
    radii_list2 = np.tile(np.array([atom.radii for atom in atom_list2]), (len(atom_list1), 1))
    radii_list1 = np.tile(np.array([[atom.radii for atom in atom_list1]]).transpose(), (1, len(atom_list2)))
    rescalar = np.add(radii_list1, radii_list2)
    dist_matrix = cdist([atom.get_coord() for atom in atom_list1], [atom.get_coord() for atom in atom_list2])
    dist_matrix = np.divide(dist_matrix, rescalar)
    return np.array(dist_matrix < threshold).ravel()


def cluster_structures(list_of_coord, cluster_method='average', clustering_number=2, metric='euclidean'):
    def fill_p_verts(p_verts, cluster):
        """
        Function: fill an empty list of lists (length equal to the number of clusters) with vertice indices according to the clustering (becomes centroids in super clustering)
        Parameters
        p_verts: empty list of lists
        cluster: cluster indices
        """

        for i, clust_idx in enumerate(cluster):
            p_idx = clust_idx - 1
            p_verts[p_idx].append(i)
        return

    zmat = hier.linkage(list_of_coord, method=cluster_method, metric=metric)
    cluster = hier.fcluster(zmat, clustering_number, criterion='maxclust', )
    clusters = [[] for _ in range(max(cluster))]
    fill_p_verts(clusters, cluster)
    #
    # plt.figure(101)
    # plt.title("ascending")
    # hier.dendrogram(zmat,color_threshold=1,p=cluster_number, truncate_mode='lastp')
    # plt.show()
    return clusters


def calculate_average_spread(clusters, distance_matrix):
    # type: (object, object) -> object
    def spread(cluster, distance_matrix):
        if len(cluster) > 1:
            return sum([distance_matrix[i][j] for i, j in combinations(cluster, 2)]) / (
                len(cluster) * (len(cluster) - 1) / 2)
        else:
            return 0

    return sum([spread(cluster, distance_matrix) for cluster in clusters]) / len(clusters)


def calculate_distance_matrix(folders, method='contact', QUITE=False):
    AS_classes = []
    for folder in folders:
        for AS_folder in [folder + '/' + f for f in os.listdir(folder) if
                          os.path.isdir(os.path.join(folder, f))]:
            if not QUITE:
                print(AS_folder)
            AS_classes.append(AlphaSpace(AS_folder))

    if method == 'contact':
        start = timeit.default_timer()

        def calculate_distance_matrix(as_class_list):
            """
            building the contact atom list for fast calculation
            """
            surface_atom_index = []
            for as_class in as_class_list:
                for atom in as_class.get_atoms(structure='surface'):
                    if str(atom.get_full_id()[3][1]) + atom.get_full_id()[4][0] not in surface_atom_index:
                        surface_atom_index.append(str(atom.get_full_id()[3][1]) + atom.get_full_id()[4][0])
            """
            building contact boolean value matrix
            """
            surface_atom_set_list = []
            ligand_atom_set_list = []
            for as_class in as_class_list:
                surface_atom_set_list.append(set([str(atom.get_full_id()[3][1]) + atom.get_full_id()[4][0] for atom in
                                                  as_class.get_atoms(structure='receptor') if
                                                  str(atom.get_full_id()[3][1]) + atom.get_full_id()[4][
                                                      0] in surface_atom_index]))
                ligand_atom_set_list.append(set([str(atom.get_full_id()[3][1]) + atom.get_full_id()[4][0] for atom in
                                                 as_class.get_atoms(structure='binder')]))
            minimal_surface_atom_set = set.intersection(*surface_atom_set_list)
            minimal_ligand_atom_set = set.intersection(*ligand_atom_set_list)
            contact_matrix_list = []
            for as_class in as_class_list:
                contact_matrix_list.append(
                    contact_detection(atom_list1=sorted([atom for atom in as_class.get_atoms(structure='receptor') if
                                                         str(atom.get_full_id()[3][1]) + atom.get_full_id()[4][
                                                             0] in minimal_surface_atom_set],
                                                        key=lambda atom: str(atom.get_full_id()[3][1]) +
                                                                         atom.get_full_id()[4][0]),
                                      atom_list2=sorted([atom for atom in as_class.get_atoms(structure='binder') if
                                                         str(atom.get_full_id()[3][1]) + atom.get_full_id()[4][
                                                             0] in minimal_ligand_atom_set],
                                                        key=lambda atom: str(atom.get_full_id()[3][1]) +
                                                                         atom.get_full_id()[4][0])))
            """
            calculate the jaccard distance matrix
            """
            return pdist(X=contact_matrix_list, metric='jaccard')

        dist_p_matrix = calculate_distance_matrix(AS_classes)
        dist_c_matrix = squareform(dist_p_matrix)
        if not QUITE:
            print("distance matrix calculated in", timeit.default_timer() - start)

    elif method == 'pro_bb_rmsd':
        imposer = Superimposer()
        protein_bb_atom_coord = [structure.get_atoms_coord(structure='receptor', backbone_only=True) for structure in
                                 AS_classes]
        dist_c_matrix = np.zeros((len(protein_bb_atom_coord), len(protein_bb_atom_coord)))
        for i in range(len(protein_bb_atom_coord)):
            for j in range(len(protein_bb_atom_coord)):
                if i < j:
                    dist_c_matrix[i][j] = imposer.set_coord(protein_bb_atom_coord[i], protein_bb_atom_coord[j])
                elif i == j:
                    dist_c_matrix[i][j] = 0.0
                else:
                    dist_c_matrix[i][j] = dist_c_matrix[j][i]
        dist_p_matrix = squareform(dist_c_matrix)

    elif method == 'face_rmsd':
        imposer = Superimposer()
        surface_atom_index = []
        for as_class in AS_classes:
            for atom in as_class.get_atoms(structure='surface'):
                if str(atom.get_full_id()[3][1]) + atom.get_full_id()[4][0] not in surface_atom_index:
                    surface_atom_index.append(str(atom.get_full_id()[3][1]) + atom.get_full_id()[4][0])
        surface_atom_index.sort()
        protein_face_atom_coord = [
            structure.get_atoms_coord(structure='receptor', backbone_only=False, atom_index=surface_atom_index) for
            structure in AS_classes]
        protein_bb_atom_coord = [structure.get_atoms_coord(structure='receptor', backbone_only=True) for structure in
                                 AS_classes]
        dist_c_matrix = np.zeros((len(protein_bb_atom_coord), len(protein_bb_atom_coord)))
        for i in range(len(protein_bb_atom_coord)):
            for j in range(len(protein_bb_atom_coord)):
                if i < j:
                    imposer.set_coord(protein_bb_atom_coord[i], protein_bb_atom_coord[j])
                    dist_c_matrix[i][j] = imposer.get_init_rms(protein_face_atom_coord[i], imposer.apply_to_coord(
                        coord_list=protein_face_atom_coord[j]))
                elif i == j:
                    dist_c_matrix[i][j] = 0.0
                else:
                    dist_c_matrix[i][j] = dist_c_matrix[j][i]
        dist_p_matrix = squareform(dist_c_matrix)

    elif method == 'face_bb_rmsd':
        imposer = Superimposer()
        surface_atom_index = []
        for as_class in AS_classes:
            for atom in as_class.get_atoms(structure='surface'):
                if str(atom.get_full_id()[3][1]) + atom.get_full_id()[4][0] not in surface_atom_index:
                    surface_atom_index.append(str(atom.get_full_id()[3][1]) + atom.get_full_id()[4][0])
        surface_atom_index.sort()
        protein_face_atom_coord = [
            structure.get_atoms_coord(structure='receptor', backbone_only=True, atom_index=surface_atom_index) for
            structure
            in AS_classes]
        protein_bb_atom_coord = [structure.get_atoms_coord(structure='receptor', backbone_only=True) for structure in
                                 AS_classes]
        dist_c_matrix = np.zeros((len(protein_bb_atom_coord), len(protein_bb_atom_coord)))
        for i in range(len(protein_bb_atom_coord)):
            for j in range(len(protein_bb_atom_coord)):
                if i < j:
                    imposer.set_coord(protein_bb_atom_coord[i], protein_bb_atom_coord[j])
                    dist_c_matrix[i][j] = imposer.get_init_rms(protein_face_atom_coord[i], imposer.apply_to_coord(
                        coord_list=protein_face_atom_coord[j]))
                elif i == j:
                    dist_c_matrix[i][j] = 0.0
                else:
                    dist_c_matrix[i][j] = dist_c_matrix[j][i]
        dist_p_matrix = squareform(dist_c_matrix)

    elif method == 'lig_rmsd':
        imposer = Superimposer()
        ligand_atom_set_list = []
        for as_class in AS_classes:
            ligand_atom_set_list.append(set([str(atom.get_full_id()[3][1]) + atom.get_full_id()[4][0] for atom in
                                             as_class.get_atoms(structure='binder')]))
        minimal_ligand_atom_set = set.intersection(*ligand_atom_set_list)
        ligand_atom_coord = [
            structure.get_atoms_coord(structure='binder', backbone_only=False, atom_index=minimal_ligand_atom_set) for
            structure in AS_classes]
        protein_bb_atom_coord = [structure.get_atoms_coord(structure='receptor', backbone_only=True) for structure in
                                 AS_classes]
        dist_c_matrix = np.zeros((len(protein_bb_atom_coord), len(protein_bb_atom_coord)))
        for i in range(len(protein_bb_atom_coord)):
            for j in range(len(protein_bb_atom_coord)):
                if i < j:
                    imposer.set_coord(protein_bb_atom_coord[i], protein_bb_atom_coord[j])
                    dist_c_matrix[i][j] = imposer.get_init_rms(ligand_atom_coord[i], imposer.apply_to_coord(
                        coord_list=ligand_atom_coord[j]))
                elif i == j:
                    dist_c_matrix[i][j] = 0.0
                else:
                    dist_c_matrix[i][j] = dist_c_matrix[j][i]
        dist_p_matrix = squareform(dist_c_matrix)

    elif method == 'lig_bb_rmsd':
        imposer = Superimposer()
        ligand_atom_set_list = []
        for i, as_class in enumerate(AS_classes):
            if not QUITE:
                print(i + 1, len(AS_classes))
            ligand_atom_set_list.append(set([str(atom.get_full_id()[3][1]) + atom.get_full_id()[4][0] for atom in
                                             as_class.get_atoms(structure='binder', backbone_only=True)]))
        minimal_ligand_atom_set = set.intersection(*ligand_atom_set_list)
        ligand_atom_coord = [
            structure.get_atoms_coord(structure='binder', backbone_only=True, atom_index=minimal_ligand_atom_set) for
            structure in AS_classes]
        protein_bb_atom_coord = [structure.get_atoms_coord(structure='receptor', backbone_only=True) for structure in
                                 AS_classes]
        dist_c_matrix = np.zeros((len(protein_bb_atom_coord), len(protein_bb_atom_coord)))
        for i in range(len(protein_bb_atom_coord)):
            for j in range(len(protein_bb_atom_coord)):
                if i < j:
                    imposer.set_coord(protein_bb_atom_coord[i], protein_bb_atom_coord[j])
                    dist_c_matrix[i][j] = imposer.get_init_rms(ligand_atom_coord[i], imposer.apply_to_coord(
                        coord_list=ligand_atom_coord[j]))
                elif i == j:
                    dist_c_matrix[i][j] = 0.0
                else:
                    dist_c_matrix[i][j] = dist_c_matrix[j][i]
                    # dist_p_matrix = squareform(dist_c_matrix)

    elif method == 'face_align':
        imposer = Superimposer()
        surface_atom_index = []
        for as_class in AS_classes:
            for atom in as_class.get_atoms(structure='surface'):
                if str(atom.get_full_id()[3][1]) + atom.get_full_id()[4][0] not in surface_atom_index:
                    surface_atom_index.append(str(atom.get_full_id()[3][1]) + atom.get_full_id()[4][0])
        """
        building contact boolean value matrix
        """
        surface_atom_set_list = []
        for as_class in AS_classes:
            surface_atom_set_list.append(set([str(atom.get_full_id()[3][1]) + atom.get_full_id()[4][0] for atom in
                                              as_class.get_atoms(structure='receptor') if
                                              str(atom.get_full_id()[3][1]) + atom.get_full_id()[4][
                                                  0] in surface_atom_index]))
        minimal_surface_atom_set = set.intersection(*surface_atom_set_list)
        protein_face_atom_coord = [
            structure.get_atoms_coord(structure='receptor', backbone_only=False, atom_index=minimal_surface_atom_set) for
            structure
            in AS_classes]
        dist_c_matrix = np.zeros((len(protein_face_atom_coord), len(protein_face_atom_coord)))
        for i in range(len(dist_c_matrix)):
            for j in range(len(dist_c_matrix)):
                if i < j:
                    dist_c_matrix[i][j] = imposer.set_coord(protein_face_atom_coord[i], protein_face_atom_coord[j])
                elif i == j:
                    dist_c_matrix[i][j] = 0.0
                else:
                    dist_c_matrix[i][j] = dist_c_matrix[j][i]
                    # dist_p_matrix = squareform(dist_c_matrix)
    return dist_c_matrix
    # From distance matrix calculate the clusters


def structure_cluster(folders, method='contact', distance_matrix=None, cluster_number=None,
                      output_average_structure=False, QUITE=False):
    """

    :param folders: list of folders where the alpha space is stored
    :param method: different methods for calculating distance matrix include: 'contact','pro_bb_rmsd','face_rmsd','face_bb_rmsd','lig_rmsd','lig_bb_rmsd'
    :param cluster_number: if you want to output the clustering result of any specific clustering number
    :return: if no cluster_number, then the penalty list of different cluster size, or the specific clustering number.
    """

    if distance_matrix is None:

        AS_classes = []
        for folder in folders:
            for AS_folder in [os.path.join(folder, f) for f in os.listdir(folder) if
                              os.path.isdir(os.path.join(folder, f))]:
                AS_classes.append(AlphaSpace(AS_folder))
        dist_c_matrix = calculate_distance_matrix(folders, method=method, QUITE=QUITE)
    else:
        dist_c_matrix = distance_matrix
    dist_p_matrix = squareform(dist_c_matrix)
    NUMBEROFSTRUCTURES = len(dist_c_matrix)

    # From distance matrix calculate the clusters
    if cluster_number is not None:
        clusters = cluster_structures(dist_p_matrix, clustering_number=cluster_number)
        if output_average_structure:
            def get_average_structure(distance_matrix, cluster):
                minimal_distance_sum = sum([distance_matrix[0][j] for j in cluster])
                minimal_distance_index = cluster[0]
                for i in cluster[1:]:
                    if sum([distance_matrix[i][j] for j in cluster]) < minimal_distance_sum:
                        minimal_distance_sum = sum([distance_matrix[i][j] for j in cluster])
                        minimal_distance_index = i
                return minimal_distance_index

            average_str_ind = [get_average_structure(dist_c_matrix, c) for c in clusters]
            return clusters, average_str_ind
        else:
            return clusters
    max_spread = calculate_average_spread(cluster_structures(dist_p_matrix, clustering_number=1),
                                          distance_matrix=dist_c_matrix)
    min_spread = calculate_average_spread(cluster_structures(dist_p_matrix, clustering_number=(NUMBEROFSTRUCTURES)),
                                          distance_matrix=dist_c_matrix)
    result = []
    for i in range((NUMBEROFSTRUCTURES)):
        average_spread = calculate_average_spread(cluster_structures(dist_p_matrix, clustering_number=i + 1),
                                                  distance_matrix=dist_c_matrix)
        penalty = (average_spread - min_spread) * ((NUMBEROFSTRUCTURES) - 2) / (max_spread - min_spread) + 2 + i
        result.append([penalty, i + 1])
    return result


def write_option_file(file_location, backupfile):
    old_opt_file = os.path.join(os.path.dirname(os.path.realpath(__file__)), backupfile)
    with open(file_location, 'w') as handle:
        handle.writelines(open(old_opt_file).readlines())


if __name__ == '__main__':
    import sys, os, pickle

    dist_mat_cache_file = None
    cwd = os.getcwd()
    sys.path.append(cwd)
    try:
        option_file_module_name = sys.argv[1][:-3]
        option = __import__(option_file_module_name)
        print("Using custom option " + option_file_module_name)
    except:
        if "option.py" not in os.listdir(cwd):
            print("Option file not found, using default options.")
            write_option_file(os.path.join(cwd, 'option.py'), "option_backup.py")
        # noinspection PyUnresolvedReferences
        import option as option
        option_file_module_name = "option"

    cld = os.path.join(cwd, option_file_module_name)
    try:
        os.mkdir(cld)
    except:
        pass


    if "dist_mat_catche.dat" in os.listdir(cld):
        print("Existing distance matrix found")
        with open(os.path.join(cld, 'dist_mat_catche.dat'), 'r') as handle:
            dist_mat_cache_file = pickle.load(handle)
    else:
        print("No distance matrix found, creating new one")
    as_master_dirs = []
    as_sub_dirs = []
    if option.alphaspace_folder_prefix is not None and len(option.alphaspace_folder_prefix) != 0:
        for folder in os.listdir(cwd):
            if folder[:len(option.alphaspace_folder_prefix)] == option.alphaspace_folder_prefix and os.path.isdir(
                    os.path.join(cwd, folder)):
                as_master_dirs.append(os.path.join(cwd, folder))
    else:
        for folder in os.listdir(cwd):
            if os.path.isdir(os.path.join(cwd, folder)):
                as_master_dirs.append(os.path.join(cwd, folder))
    for as_master_dir in as_master_dirs:
        for as_sub_dir in os.listdir(as_master_dir):
            if os.path.isdir(os.path.join(as_master_dir, as_sub_dir)):
                as_sub_dirs.append(os.path.join(as_master_dir, as_sub_dir))

    # method = 'lig_bb_rmsd'
    # folders = ['./MD/Design/as1','./MD/Design/as2','./MD/Design/as3','./MD/Design/as4','./MD/Design/as5']
    # directories = []
    # for folder in folders:
    #     directories += ([os.path.join(folder, f) for f in os.listdir(folder) if
    #                      os.path.isdir(os.path.join(folder, f))])
    # with open('./MD/Design/' + method, 'rb') as handle:
    #     distance_matrix = pickle.load(handle)
    if dist_mat_cache_file is None:
        dist_mat_cache_file = calculate_distance_matrix(folders=as_master_dirs, method=option.clustering_method,
                                                        QUITE=option.QUITE)
        with open(os.path.join(cld, 'dist_mat_catche.dat'), 'w') as handle:
            pickle.dump(dist_mat_cache_file, handle)
    clusters, cluster_rep = structure_cluster(as_master_dirs, distance_matrix=dist_mat_cache_file,
                                              method=option.clustering_method, output_average_structure=True,
                                              cluster_number=option.cluster_number, QUITE=True)

    with open(os.path.join(cld, 'cluster_result.txt'), 'w') as handle:
        for i in range(option.cluster_number):
            handle.write("Cluster {} ({} structures):\n".format(i, len(clusters[i])))
            print(clusters[i], cluster_rep[i])
            for c_ind in clusters[i]:
                if c_ind == cluster_rep[i]:
                    handle.write("*R* {}\n".format(as_sub_dirs[c_ind]))
                else:
                    handle.write("    {}\n".format(as_sub_dirs[c_ind]))
    with open(os.path.join(cld, 'cluster_time_ev.txt'), 'w') as handle:
        ts_list = [0 for _ in range(len(as_sub_dirs))]
        for cluster_ind in range(len(clusters)):
            for i in clusters[cluster_ind]:
                ts_list[i] = cluster_ind
        for i in range(len(ts_list)):
            handle.write("{}\t{}\n".format(i, ts_list[i]))
