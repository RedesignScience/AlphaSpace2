from itertools import combinations_with_replacement

import numpy as np
from mdtraj.geometry import _geometry
from mdtraj.geometry.sasa import _ATOMIC_RADII
from mdtraj.core import element
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial import Voronoi, Delaunay, cKDTree
from scipy.spatial.distance import cdist
from scipy.spatial.distance import squareform


def _getTetrahedronVolume(coord_list: list):
    """
    Calculate the volume of a tetrahedron described by four 3-d points.
    :param coord_list: list 4*3 4 points xyz
    :return: float
    """
    coord_matrix = np.concatenate((np.array(coord_list), np.ones((4, 1))), axis=1)
    volume = np.abs(np.linalg.det(coord_matrix) / 6)

    return volume


def _getTetrahedronVolumes(coord_group_list):
    volumes = np.empty((len(coord_group_list)), dtype=coord_group_list.dtype)

    for i in range(len(coord_group_list)):
        volumes[i] = _getTetrahedronVolume(coord_group_list[i])

    return volumes

def _getNonpolarRatio(protein_snapshot, alpha_lining_idx):
    areas = getSASA(protein_snapshot=protein_snapshot)


    is_nonpolar = np.array(
        [atom.element not in [element.oxygen,element.nitrogen,element.sulfur]
         for atom in protein_snapshot.topology.atoms]).astype(int)

    nonpolar_ratio = (areas * is_nonpolar)[alpha_lining_idx].sum(axis=1).astype(float) / (
            areas[alpha_lining_idx].astype(float).sum(axis=1) + 0.00001)

    return nonpolar_ratio


def _getContactMatrix(coord_list_1, coord_list_2, threshold):
    """
    For two sets of points A and B, generate the contact matrix M,
    where M(i,j) = (if Ai,Bj is in contact)
    get M by N bool matrix of if there is a contact.
    :param coord_list_1: np.ndarray N * 3
    :param coord_list_2: np.ndarray M * 3
    :param threshold: float
    :return: np.ndarray N * M
    """
    distance_matrix = cdist(coord_list_1, coord_list_2)
    return (distance_matrix < threshold).astype(int)


def _getIfContact(checked_coord_list, ref_coord_list, threshold):
    """
    Check which one in the coordinate list is in contact with the second coordinate
    :param checked_coord_list: list of array N*3 to be checked
    :param ref_coord_list: list of array M*3
    :param threshold: float
    :return: np.ndarray 2*N [index in contact][index not in contact]
    """
    return np.where(_getContactMatrix(checked_coord_list, ref_coord_list, threshold))


def _getGridVolume(coord_list, threshold=1.6, resolution=0.05):
    """
    Calculate the volume of a point set using grid point approximation
    :param coord_list: array N * 3
    :param threshold: float
    :param resolution: float
    :return: float
    """
    max_coord = np.max(coord_list, axis=0)
    min_coord = np.min(coord_list, axis=0)
    coord_range = np.array([min_coord, max_coord]).transpose().tolist()
    x, y, z = [np.arange(start=ax[0] - threshold, stop=ax[1] + threshold, step=resolution) for ax in coord_range]
    grid_coords = np.array(np.meshgrid(x, y, z)).transpose().reshape((-1, 3))

    grid_count = len(_getIfContact(grid_coords, coord_list, threshold=threshold)[0])
    return grid_count * (resolution ** 3)


def _getCos(v1, v2):
    """
    Calculate the Cos of the angle be vector v1 and v2
    :param v1: numpy.ndarray vector1
    :param v2: numpy.ndarray vector2
    :return: float
    """

    def unit_vector(vector):
        norm = np.linalg.norm(vector)
        assert norm > 0
        """ Returns the unit vector of the vector.  """
        return vector / norm

    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.clip(np.dot(v1_u, v2_u), -1.0, 1.0)


def combination_intersection_count(indices_list: list, total_index: int) -> np.ndarray:
    """

    Given a list of indices list, such as [ [1,2,3], [4,5,6] , [2,3,4]]
    This function calculates the intersection count between each pair in the list
    if indices_list is a N list of M indices, the return array D is a N * N ndarray, where Dij is the intersection count
    between i and j.
    Note this D matrix is symmetrical.
    The total_index is maximum of all indices in the indices list.
    :param indices_list: list [[indices].]
    :param total_index: int, upper limits of the indices
    :return: np.ndarray N * N where (i,j) means count of i and j element wise intersection


    Parameters
    ----------
    indices_list : list [[indices].]
    total_index : total_index: int, upper limits of the indices

    Returns
    -------
    intersection binary matrix : np.ndarray
         N * N where (i,j) means count of i and j element wise intersection

    """
    item_binary_vectors = np.empty((len(indices_list), total_index))
    item_binary_vectors.fill(0)
    for pocket_idx, lining_atoms_idx in enumerate(indices_list):
        item_binary_vectors[pocket_idx].put(lining_atoms_idx, 1)

    overlap_matrix = np.empty((item_binary_vectors.shape[0], item_binary_vectors.shape[0]))
    overlap_matrix.fill(0)
    count = 0
    for i, j in combinations_with_replacement(range(len(indices_list)), 2):
        print(count, len(overlap_matrix) ** 2 / 2)
        overlap = np.dot(item_binary_vectors[i], item_binary_vectors[j])
        overlap_matrix[i][j] = overlap_matrix[j][i] = overlap

    return overlap_matrix


def combination_union_count(indices_list, total_index: int) -> np.ndarray:
    """
    Given a list of indices list, such as [ [1,2,3], [4,5,6] , [2,3,4]]
    This function calculates the union count between each pair in the list
    if indices_list is a N list of M indices, the return array D is a N * N ndarray, where Dij is the union count
    between i and j.
    Note this D matrix is symmetrical.
    The total_index is maximum of all indices in the indices list.
    :param indices_list: list [[indices].]
    :param total_index: int, upper limits of the indices
    :return: np.ndarray N * N where (i,j) means count of i and j element wise union
    """

    indices_list = np.array(indices_list)
    item_binary_vectors = np.empty((len(indices_list), total_index))
    item_binary_vectors.fill(0)
    for pocket_idx, lining_atoms_idx in enumerate(indices_list):
        item_binary_vectors[pocket_idx].put(lining_atoms_idx, 1)

    intersection_matrix = np.empty((item_binary_vectors.shape[0], item_binary_vectors.shape[0]))
    intersection_matrix.fill(0)
    for i, j in combinations_with_replacement(range(len(indices_list)), 2):
        intersection = np.count_nonzero(item_binary_vectors[i] + item_binary_vectors[j])
        intersection_matrix[i][j] = intersection_matrix[j][i] = intersection

    return intersection_matrix


def count_intersect(a, b):
    return np.count_nonzero(a + b)



def getSASA(protein_snapshot, cover_atom_coords=None):
    """
    Calculate the absolute solvent accessible surface area.
    First calculate the SASA of the receptor by itself, then subtract it with sasa with the AAC.
    AAC are set to resemble Carbon with a radii - 0.17

    """
    probe_radius = 0.14
    n_sphere_points = 960

    if cover_atom_coords is None:
        xyz = np.array(protein_snapshot.xyz, dtype=np.float32)
        atom_radii = [_ATOMIC_RADII[atom.element.symbol] for atom in protein_snapshot.topology.atoms]
    else:
        xyz = np.array(np.expand_dims(np.concatenate((protein_snapshot.xyz[0], cover_atom_coords), axis=0), axis=0),
                       dtype=np.float32)
        atom_radii = [_ATOMIC_RADII[atom.element.symbol] for atom in protein_snapshot.topology.atoms] + [0.17 for _ in
                                                                                                         range(
                                                                                                             xyz.shape[
                                                                                                                 1])]
    radii = np.array(atom_radii, np.float32) + probe_radius
    atom_mapping = np.arange(xyz.shape[1], dtype=np.int32)
    out = np.zeros((1, xyz.shape[1]), dtype=np.float32)
    _geometry._sasa(xyz, radii, int(n_sphere_points), atom_mapping, out)
    return out[:, :protein_snapshot.xyz.shape[1]][0]


def screenContact(data, binder_xyz, threshold):
    """
    Mark the contact in AS_Data as true for each frame.
    :param data: Snapshot
    :param binder_xyz: np.ndarray
    :param threshold: float
    """
    assert len(binder_xyz.shape) == 3
    assert binder_xyz.shape[-1] == 3

    alpha_xyz = data.xyz()
    contact_matrix = _getContactMatrix(alpha_xyz, binder_xyz[data.snapshot_idx()], threshold)
    data[:, 12] = contact_matrix.any(axis=1).astype(int)


def _tessellation(**kwargs):
    """
    This is the main AlphaSpace function, it's self contained so you can run it in
    multiprocessing module.


    """

    receptor_xyz = kwargs['receptor_xyz']
    binder_xyz = kwargs['binder_xyz']
    atom_radii = kwargs['atom_radii']

    config = kwargs['config']
    snapshot_idx = kwargs['snapshot_idx']

    try:
        cluster_method = config.cluster_method
    except:
        cluster_method = 'average_linkage'

    # Generate Raw Tessellation simplexes
    raw_alpha_lining_idx = Delaunay(receptor_xyz).simplices
    # Take coordinates from xyz file
    raw_alpha_lining_xyz = np.take(receptor_xyz, raw_alpha_lining_idx[:, 0].flatten(), axis=0)

    # generate alpha atom coordinates
    raw_alpha_xyz = Voronoi(receptor_xyz).vertices
    # Calculate alpha sphere radii
    raw_alpha_sphere_radii = np.linalg.norm(raw_alpha_lining_xyz - raw_alpha_xyz, axis=1)

    # Filter the data based on radii cutoff
    filtered_alpha_idx = np.where(np.logical_and(config.min_r<= raw_alpha_sphere_radii,
                                                 raw_alpha_sphere_radii <= config.max_r))[0]

    filtered_alpha_radii = np.take(raw_alpha_sphere_radii, filtered_alpha_idx)

    alpha_lining = np.take(raw_alpha_lining_idx, filtered_alpha_idx, axis=0)

    filtered_alpha_xyz = np.take(raw_alpha_xyz, filtered_alpha_idx, axis=0)

    if cluster_method == 'average_linkage':
        # cluster the remaining vertices to assign index of belonging pockets
        zmat = linkage(filtered_alpha_xyz, method='average')

        alpha_pocket_index = fcluster(zmat, config.clust_dist / 10,
                                      criterion='distance') - 1  # because cluster index start from 1
    elif cluster_method == 'hdbscan':
        import hdbscan
        clusterer = hdbscan.HDBSCAN(metric='euclidean', min_samples=config.hdbscan_min_samples)
        clusterer.fit(filtered_alpha_xyz)
        alpha_pocket_index = clusterer.labels_

    else:
        raise Exception('Known Clustering Method: {}'.format(cluster_method))

    # Load trajectories
    filtered_lining_xyz = np.take(receptor_xyz, alpha_lining, axis=0)
    # calculate the polarity of alpha atoms
    _total_space = np.array(
        [_getTetrahedronVolume(i) for i in filtered_lining_xyz]) * 1000  # here the 1000 is to convert nm^3 to A^3

    _nonpolar_space = _polar_space = _total_space / 2

    if binder_xyz is not None:

        """
        Calculate the contact matrix, and link each alpha with closest atom
        """
        dist_matrix = cdist(filtered_alpha_xyz, binder_xyz)

        min_idx = np.argmin(dist_matrix, axis=1)
        mins = np.min(dist_matrix, axis=1) * 10  # nm to A
        is_contact = mins < config.hit_dist

    else:
        min_idx = np.zeros(filtered_alpha_xyz.shape[0])
        mins = np.zeros(filtered_alpha_xyz.shape[0])
        is_contact = np.zeros(filtered_alpha_xyz.shape[0])

    """lining atom asa"""
    _xyz = np.array(np.expand_dims(receptor_xyz, axis=0),
                    dtype=np.float32)
    dim1 = _xyz.shape[1]
    atom_mapping = np.arange(dim1, dtype=np.int32)
    asa = np.zeros((1, dim1), dtype=np.float32)
    radii = np.array(atom_radii, np.float32) + config.probe_radius
    _geometry._sasa(_xyz, radii, int(config.n_sphere_points), atom_mapping, asa)

    alpha_lining_asa = np.take(asa[0], alpha_lining).sum(axis=1) * 100  # nm2 to A2

    """set contact to active if use ligand contact is True"""
    is_active = is_contact if config.screen_by_lig_cntct else np.zeros_like(alpha_pocket_index)

    data = np.concatenate((np.array([range(alpha_pocket_index.shape[0])]).transpose(),  # 0         idx
                           np.full((alpha_pocket_index.shape[0], 1), snapshot_idx, dtype=int),  # 1         snapshot_idx
                           filtered_alpha_xyz,  # 2 3 4     x y z
                           alpha_lining,  # 5 6 7 8   lining_atom_idx_1 - 4
                           np.expand_dims(_polar_space, axis=1),  # 9         polar_space 0
                           np.expand_dims(_nonpolar_space, axis=1),  # 10        nonpolar_space 0
                           np.expand_dims(is_active, axis=1),  # 11        is_active 1
                           np.expand_dims(is_contact, axis=1),  # 12        isContact 0
                           np.expand_dims(alpha_pocket_index, axis=1),  # 13        pocket_idx
                           np.expand_dims(filtered_alpha_radii, axis=1),  # 14        radii
                           np.expand_dims(min_idx, axis=1),  # 15        closest atom idx
                           np.expand_dims(mins, axis=1),  # 16        closest atom dist
                           np.expand_dims(alpha_lining_asa, axis=1)  # 17 total lining atom asa
                           ), axis=-1)

    print('{} snapshot processed'.format(snapshot_idx + 1))
    return data

def cluster_by_overlap(vectors, total_index, overlap_cutoff):
    """
    Cluster a list of binary vectors based on lining atom overlap

    Parameters
    ----------
    vectors
    total_index
    overlap_cutoff

    Returns
    -------

    """
    # calculate jaccard_diff_matrix
    intersection_matrix = combination_intersection_count(vectors, total_index)

    union_matrix = combination_union_count(vectors, total_index)

    jaccard_diff_matrix = 1 - intersection_matrix / union_matrix

    cluster_index = list(fcluster(Z=linkage(squareform(jaccard_diff_matrix), method='average'),
                                  t=overlap_cutoff,
                                  criterion='distance') - 1)

    cluster_list = {i: [] for i in range(max(cluster_index) + 1)}

    for cluster_i, item_idx in enumerate(cluster_index):
        cluster_list[item_idx].append(cluster_i)

    return cluster_list


def _bestProbeType(beta_atom):
    """
    Parameters
    ----------
    beta_atom : AS_BetaAtom

    Get the probe type for the best score in this beta atom.

    Returns
    -------
    probe_type : str
            ['C', 'Br', 'F', 'Cl', 'I', 'OA', 'SA', 'N', 'P']
    """
    _best_score_index = min(range(9), key=lambda i: beta_atom.score)

    return ['C', 'Br', 'F', 'Cl', 'I', 'OA', 'SA', 'N', 'P'][_best_score_index]


def _isPocketConnected(p1, p2):
    """
    Check if two pockets are connected

    Parameters
    ----------
    p1 : AS_Pocket
    p2 : AS_Pocket

    Returns
    -------
    bool

    """
    if set(p1.lining_atoms_idx).intersection(p2.lining_atoms_idx):
        pocket_vector1 = p1.lining_atoms_centroid - p1.centroid
        pocket_vector2 = p2.lining_atoms_centroid - p2.centroid
        if _getCos(pocket_vector1, pocket_vector2) > 0:  # pocket vector facing inwards
            return True
    return False


def _prune_dpockets(d_pocket_dict, sample_ratio=1.0):
    """

    This is used in pruning pockets in dpockets and selecting a subset for leader follower clustering.

    Parameters
    ----------
    d_pocket_dict : dict
    sample_ratio: float

    Returns
    -------

    """

    leader = []
    labels = []
    assert sample_ratio <= 1.0
    for d_pocket_idx, pockets in d_pocket_dict.items():
        n_leaders = int(np.ceil(len(pockets) * float(sample_ratio)))
        leader.extend(np.random.choice(pockets, n_leaders, replace=False))
        labels.extend([d_pocket_idx] * n_leaders)
    return leader, labels


def _group(label_list):
    d = [[] for _ in range(max(label_list) + 1)]
    for i, item in enumerate(label_list):
        d[item].append(i)
    return d


def _findInRange(query_points, ref_points, cutoff):
    indices = np.where(cdist(query_points, ref_points) <= cutoff)[0]
    indices = np.unique(indices)
    return indices


def _markInRange(query_points, ref_points, cutoff):
    indices = _findInRange(query_points, ref_points, cutoff)
    query_bool = np.zeros(len(query_points), dtype=np.bool)
    query_bool[indices] = 1
    return query_bool


def _binCluster(coords, bin_size, bin_buffer_ratio=0.01, distance=5.4):
    """
    Bin the coordinates first and cluster based on bin centroids

    Parameters
    ----------
    coords:
            np.array shape = (n,3)

    bin_size: [x_bin_size,y_bin_size,z_bin_size]
            list

    bin_buffer_ratio: float

    Returns
    -------

    labels : label for all points
            np.array shape = (n)

    """
    # find bin border with buffer ratio
    xyz_max = np.max(coords, 0)
    xyz_min = np.min(coords, 0)

    buffer_size = (xyz_max - xyz_min) * bin_buffer_ratio

    xyz_max += buffer_size
    xyz_min -= buffer_size

    # create bin arrays
    xyz_bins = []
    binned_idx = []
    bin_number = [None, None, None]
    for dim in range(3):
        bin_number[dim] = int((xyz_max[dim] - xyz_min[dim]) / float(bin_size[dim])) + 1

        xyz_bins.append(np.linspace(xyz_min[dim], xyz_max[dim], bin_number[dim] + 1))
        val = coords[:, dim]
        ind = np.digitize(val, xyz_bins[dim]) - 1

        binned_idx.append(ind)

    binned_idx = np.array(binned_idx)

    add_location = np.transpose(binned_idx)

    cube_coord_sum = np.zeros(bin_number + [3])

    cube_idx_count = np.zeros(bin_number)

    for idx in np.arange(len(coords)):
        # print(cube_coord_sum[add_location[idx]])
        add_idx = add_location[idx]
        cube_coord_sum[add_idx[0], add_idx[1], add_idx[2]] += coords[idx]
        cube_idx_count[add_idx[0], add_idx[1], add_idx[2]] += 1

    non_zero_cube_idx = np.where(cube_idx_count > 0)
    non_zero_cube_coord = np.divide(cube_coord_sum[non_zero_cube_idx],
                                    np.expand_dims(cube_idx_count[non_zero_cube_idx], 1))
    zmat = linkage(non_zero_cube_coord, method='average')
    label = fcluster(zmat, distance, criterion='distance') - 1

    tree = cKDTree(non_zero_cube_coord)
    dist, ind = tree.query(coords, 1, n_jobs=1)

    pocket_label = label[ind]

    return pocket_label
