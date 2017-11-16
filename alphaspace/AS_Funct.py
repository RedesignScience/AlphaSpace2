import numpy as np
from mdtraj.geometry import _geometry
from mdtraj.geometry.sasa import _ATOMIC_RADII
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial import Voronoi, Delaunay
from scipy.spatial.distance import cdist


def getTetrahedronVolume(coord_list):
    """
    generate the volume
    :param coord_list: list 4*3
    :return: float
    """
    coord_matrix = np.concatenate((np.array(coord_list), np.ones((4, 1))), axis=1)
    volume = np.abs(np.linalg.det(coord_matrix) / 6)

    return volume


def getContactMatrix(coord_list_1, coord_list_2, threshold):
    """
    get M by N bool matrix of if there is a contact.
    :param coord_list_1: np.ndarray N * 3
    :param coord_list_2: np.ndarray M * 3
    :param threshold: float
    :return: np.ndarray N * M
    """
    distance_matrix = cdist(coord_list_1, coord_list_2)
    return (distance_matrix < threshold).astype(int)


def getIfContact(coord_list_1, coord_list_2, threshold):
    """
    check which one in the coordinate list is in contact with the second coordinate
    :param coord_list_1: list of array N*3
    :param coord_list_2: list of array M*3
    :param threshold: float
    :return: np.ndarray 2*N
    """

    return np.where(getContactMatrix(coord_list_1, coord_list_2, threshold))


def getGridVolume(coord_list, threshold=1.6, resolution=0.05):
    """
    calculate the volume of a point set using grid point approximation
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

    grid_count = len(getIfContact(grid_coords, coord_list, threshold=threshold)[0])
    return grid_count * (resolution ** 3)


def getSASA(protein_snapshot, cover_atom_coords=None):
    """
    Calculate the absolute solvent accessible surface area.
    First calculate the SASA of the receptor by itself, then subtract it with sasa with the AAC.
    AAC are set to resemble Carbon with a radii - 0.17
    :param protein_snapshot: mdtraj object
    :param cover_atom_coords: np.ndarray n*3
    :return: np.ndarray n
    """
    probe_radius = 0.14
    n_sphere_points = 960

    if cover_atom_coords is None:
        xyz = np.array(protein_snapshot.xyz,dtype=np.float32)
        atom_radii = [_ATOMIC_RADII[atom.element.symbol] for atom in protein_snapshot.topology.atoms]
        # atom_mapping = np.arange(xyz.shape[1],dtype=np.int32)
        # out = np.zeros((1,xyz.shape[1]),dtype=np.float32)
        # radii = np.array(atom_radii,np.float32) + probe_radius
        # _geometry._sasa(xyz,radii,int(n_sphere_points),atom_mapping,out)
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
    out = out[:, :protein_snapshot.xyz.shape[1]]
    return out[0]


def screenContact(data, binder_xyz, threshold):
    """
    Mark the contact in AS_Data as true for each frame.
    :param data: AS_Data
    :param binder_xyz: np.ndarray
    :param threshold: float
    """
    assert len(binder_xyz.shape) == 3
    assert binder_xyz.shape[-1] == 3

    for snapshot_idx in data.snapshots_idx:
        alpha_idx = data.snapshot_alpha_idx(snapshot_idx)
        alpha_xyz = data.xyz(alpha_idx)
        contact_matrix = getContactMatrix(alpha_xyz, binder_xyz[snapshot_idx], threshold)
        data[alpha_idx, 12] = contact_matrix.any(axis=1).astype(int)



def _tessellation(queue,arglist):
    assert len(arglist) == 4

    protein_snapshot, config, snapshot_idx, is_polar = arglist
    xyz = protein_snapshot.xyz[0]
    # Generate Raw Tessellation simplexes
    raw_alpha_lining_idx = Delaunay(xyz).simplices
    # Take coordinates from xyz file
    raw_alpha_lining_xyz = np.take(xyz, raw_alpha_lining_idx[:, 0].flatten(), axis=0)

    # generate alpha atom coordinates
    raw_alpha_xyz = Voronoi(xyz).vertices
    # Calculate alpha sphere radii
    raw_alpha_sphere_radii = np.linalg.norm(raw_alpha_lining_xyz - raw_alpha_xyz,axis=1)

    # Filter the data based on radii cutoff
    filtered_alpha_idx = np.where(np.logical_and(config.min_r / 10.0 <= raw_alpha_sphere_radii,
                                                 raw_alpha_sphere_radii <= config.max_r / 10.0))[0]

    filtered_alpha_radii = np.take(raw_alpha_sphere_radii,filtered_alpha_idx)

    alpha_lining = np.take(raw_alpha_lining_idx,filtered_alpha_idx,axis=0)

    filtered_alpha_xyz = np.take(raw_alpha_xyz,filtered_alpha_idx,axis=0)

    # cluster the remaining vertices to assign index of belonging pockets
    zmat = linkage(filtered_alpha_xyz,method='average')

    alpha_pocket_index = fcluster(zmat, config.clust_dist / 10,
                                  criterion='distance') - 1  # because cluster index start from 1

    # Load trajectories
    filtered_lining_xyz = np.take(xyz, alpha_lining, axis=0)
    # calculate the polarity of alpha atoms
    _total_score = np.array(
            [getTetrahedronVolume(i) for i in
             filtered_lining_xyz]) * 1000  # here the 1000 is to convert nm^3 to A^3

    element = [str(atom.element.symbol) for atom in protein_snapshot.topology._atoms]
    atom_radii = [_ATOMIC_RADII[e] for e in element]
    alpha_radii = [0.17 for _ in range(len(alpha_pocket_index))]

    """"""

    _xyz = np.array(np.expand_dims(np.concatenate((xyz, filtered_alpha_xyz), axis=0), axis=0),
                    dtype=np.float32)
    dim1 = _xyz.shape[1]
    atom_mapping = np.arange(dim1, dtype=np.int32)
    covered = np.zeros((1, dim1), dtype=np.float32)
    radii = np.array(atom_radii + alpha_radii, np.float32) + config.probe_radius
    _geometry._sasa(_xyz, radii, int(config.n_sphere_points), atom_mapping, covered)
    covered = covered[:, :xyz.shape[0]]

    _xyz = np.array(np.expand_dims(xyz, axis=0),
                    dtype=np.float32)
    dim1 = _xyz.shape[1]
    atom_mapping = np.arange(dim1, dtype=np.int32)
    total = np.zeros((1, dim1), dtype=np.float32)
    radii = np.array(atom_radii, np.float32) + config.probe_radius
    _geometry._sasa(_xyz, radii, int(config.n_sphere_points), atom_mapping, total)

    atom_sasa = (total - covered)[0] * 100  # nm^2 to A^2 convertion

    """"""
    # atom_sasa = get_abs_sasa(xyz,filtered_alpha_xyz,np.array(atom_radii+alpha_radii,dtype=np.float32))
    pocket_sasa = np.take(atom_sasa,alpha_lining)

    polar_ratio = np.average(np.take(is_polar.astype(float),alpha_lining),axis=1,weights=pocket_sasa)

    _polar_score = _total_score * polar_ratio

    _nonpolar_score = _total_score - _polar_score

    """
    Calculate the contact matrix, and link each alpha with closest atom.
    """

    dist_matrix = cdist(filtered_alpha_xyz, xyz)
    min_idx = np.argmin(dist_matrix, axis=1)
    mins = np.min(dist_matrix, axis=1) * 10  # nm to A
    is_contact = mins < config.hit_dist

    is_active = is_contact if config.screen_by_lig_cntct else np.ones_like(alpha_pocket_index)

    data = np.concatenate((np.zeros((len(alpha_pocket_index), 1)),  # 0         idx
                           np.full((len(alpha_pocket_index), 1), snapshot_idx),  # 1         snapshot_idx
                           filtered_alpha_xyz,  # 2 3 4     x y z
                           alpha_lining,  # 5 6 7 8   lining_atom_idx_1 - 4
                           np.expand_dims(_polar_score, axis=1),  # 9         polar_score 0
                           np.expand_dims(_nonpolar_score, axis=1),  # 10        nonpolar_score 0
                           np.expand_dims(is_active, axis=1),  # 11        is_active 1
                           np.expand_dims(is_contact, axis=1),  # 12        is_contact 0
                           np.expand_dims(alpha_pocket_index, axis=1),  # 13        pocket_idx
                           np.expand_dims(filtered_alpha_radii, axis=1),  # 14        radii
                           np.expand_dims(min_idx, axis=1),  # 15        closest atom idx
                           np.expand_dims(mins, axis=1),  # 16        closest atom dist
                           ),axis=-1)

    assert data.shape[1] == 17

    """
    0       idx
    1       snapshot_idx
    2       x
    3       y
    4       z
    5       lining_atom_idx_1
    6       lining_atom_idx_1
    7       lining_atom_idx_1
    8       lining_atom_idx_1
    9       polar_score 0
    10      nonpolar_score 0
    11      is_active 1
    12      is_contact 0
    13      pocket_idx
    14      radii
    15      closest atom idx
    16      closest atom dist
    """
    print('{} snapshot processed'.format(snapshot_idx))
    queue.put(data)
    return
