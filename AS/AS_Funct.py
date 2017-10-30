import numpy as np
from scipy.spatial.distance import cdist
from AS_Config import AS_Config
from scipy.spatial import Voronoi, Delaunay
from scipy.cluster.hierarchy import linkage, fcluster
from mdtraj.core.topology import Residue,Atom
from mdtraj.core.element import *


def _tessellation(cluster):
    """
    perform tessellation in order to generate the cluster of alpha atoms.
    :param cluster: object
    """

    # Generate Raw Tessellation simplices
    raw_alpha_lining_idx = Delaunay(cluster.receptor_snapshot.xyz[0]).simplices
    # Take coordinates from xyz file
    raw_alpha_lining_xyz = np.take(cluster.receptor_snapshot.xyz[0], raw_alpha_lining_idx[:, 0].flatten(), axis=0)

    # generate alpha atom coordinates
    raw_alpha_xyz = Voronoi(cluster.receptor_snapshot.xyz[0]).vertices

    # Calculate alpha sphere radii
    raw_alpha_sphere_radii = np.linalg.norm(raw_alpha_lining_xyz - raw_alpha_xyz, axis=1)

    # Filter the data based on radii cutoff
    filtered_alpha_idx = np.where(np.logical_and(cluster.config.min_r / 10.0 <= raw_alpha_sphere_radii,
                                                 raw_alpha_sphere_radii <= cluster.config.max_r / 10.0))[0]

    cluster.alpha_lining = np.take(raw_alpha_lining_idx, filtered_alpha_idx, axis=0)

    filtered_alpha_xyz = np.take(raw_alpha_xyz, filtered_alpha_idx, axis=0)

    # cluster the remaining vertices to assign index of belonging pockets
    zmat = linkage(filtered_alpha_xyz, method='average')
    pocket_idx = fcluster(zmat, cluster.config.pocket_cluster_distance / 10, criterion='distance')  # /10 turn A to nm

    cluster.alpha_pocket_index = pocket_idx - 1  # because cluster index start from 1

    # print(max(pocket_idx),min(pocket_idx))

    # Reorganize into list of pockets
    # self.pocket_alpha_atoms = [[] for _ in range(max(cluster))]
    # for alpha_cluster_i, alpha_atom_idx in sorted(
    #         zip(self.alpha_pocket_index, range(len(self.alpha_pocket_index)))):
    #     self.pocket_alpha_atoms[alpha_cluster_i].append(alpha_atom_idx)

    # Generate Residue container for pockets and add in atoms as AAC
    for _ in range(max(pocket_idx) - 1):
        residue = cluster.top.add_residue(name='ASC', chain=cluster.top.chain(0))
        residue.cluster = cluster

    for i, pocket_index in enumerate(cluster.alpha_pocket_index):
        atom = cluster.top.add_atom('AAC', Element(), cluster.top.residue(pocket_index), pocket_index)
        atom.index = i

    update_atom_methods(Atom)
    update_residue_method(Residue)

    for pocket in cluster.pockets:
        pocket.cluster = cluster  # assign the parent cluster
        alpha_index = [alpha.index for alpha in pocket.atoms]
        pocket.lining_atom_idx = cluster._get_lining_atoms(alpha_index)  # Assign pocket lining atoms
        pocket.lining_residue_idx = cluster._get_lining_residues(alpha_index)
    # Load trajectories
    cluster.traj.xyz = np.expand_dims(filtered_alpha_xyz, axis=0)
    filtered_lining_xyz = np.take(cluster.receptor_snapshot.xyz[0], cluster.alpha_lining, axis=0)
    # calculate the polarity of alpha atoms
    cluster._total_score = np.array([getTetrahedronVolume(i) for i in filtered_lining_xyz])

    pocket_sasa = np.take(cluster._get_SASA(), cluster.alpha_lining)

    polar_ratio = np.average(np.take(cluster.is_polar, cluster.alpha_lining), axis=1, weights=pocket_sasa)

    cluster._polar_score = cluster._total_score * polar_ratio

    cluster._nonpolar_score = cluster._total_score - cluster._polar_score

    cluster._contact = np.zeros(cluster.top.n_atoms, dtype=int)

    cluster.is_run = True
    return cluster

def getTetrahedronVolume(coord_list):
    """
    generate the volume
    :param coord_list: list N*4*3
    :return: float
    """
    coord_matrix = np.concatenate((np.array(coord_list), np.ones((4, 1))), axis=1)

    return np.abs(np.linalg.det(coord_matrix) / 6)


def getContactMatrix(coord_list_1, coord_list_2, threshold=None):
    if threshold is None:
        threshold = AS_Config().contact_threshold
    distance_matrix = cdist(coord_list_1, coord_list_2)
    return distance_matrix < threshold


def checkContact(coord_list_1, coord_list_2, threshold):
    """
    check which one in the coordinate list is in contact with the second coordinate
    :param coord_list_1: list of array N*3
    :param coord_list_2: list of array M*3
    :param threshold: float
    :return: np.array 2*N
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

    grid_count = len(checkContact(grid_coords, coord_list, threshold=threshold)[0])
    return grid_count * (resolution ** 3)


def update_atom_methods(target):
    def get_pocket(self):
        return self.residue

    def set_cluster(self, cluster):
        self.cluster = cluster

    @property
    def cluster(self):
        return self.residue.cluster

    def get_polar_score(self):
        return self.cluster._polar_score[self.index]

    def get_nonpolar_score(self):
        return self.cluster._nonpolar_score[self.index]

    def get_total_score(self):
        return self.cluster._total_score[self.index]

    def get_contact(self):
        return self.contact

    def set_contact(self, contact):
        self.contact = contact



    target.get_pocket = get_pocket
    target.set_cluster = set_cluster
    target.cluster = cluster
    target.get_polar_score = get_polar_score
    target.get_nonpolar_score = get_nonpolar_score
    target.get_contact = get_contact
    target.set_contact = set_contact


def update_residue_method(target):
    @property
    def alphas(self):
        return self.atoms
    target.alphas = alphas

    def get_polar_score(self):
        return np.sum([alpha.get_polar_score() for alpha in self.alphas])
    target.get_polar_score = get_polar_score


    def get_nonpolar_score(self):
        return np.sum([alpha.get_nonpolar_score() for alpha in self.alphas])
    target.get_nonpolar_score = get_nonpolar_score

    def get_total_score(self):
        return np.sum([alpha.get_total_score() for alpha in self.alphas])
    target.get_total_score = get_total_score


    @property
    def parent(self):
        return self.cluster
    target.parent = parent


    def get_contact(self):
        return np.any([alpha.get_contact() for alpha in self.alphas])
    target.get_contact = get_contact


    def get_alpha_index(self):
        return [atom.index for atom in self.alphas]
    target.get_alpha_index = get_alpha_index

    def get_lining_atoms(self):
        return self.cluster._get_pocket_lining_atoms(self)
        # return self.cluster._get_pocket_lining_atoms(self)
    target.get_lining_atoms = get_lining_atoms

    def get_lining_residues(self):
        return self.cluster._get_pocket_lining_residues(self)
    target.get_lining_residues = get_lining_residues

    def get_centroid(self):
        return self.cluster._get_cluster_centroid(list(self.get_alpha_index))
    target.get_centroid = get_centroid



