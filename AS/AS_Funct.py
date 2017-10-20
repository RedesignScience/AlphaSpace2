import numpy as np
from scipy.spatial.distance import cdist
from AS.AS_Config import AS_Config


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
    print(coord_list_1.shape, coord_list_2.shape)
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

    def get_cluster(self):
        return self.cluster

    def get_polar_score(self):
        return self.cluster.polar_score[self.index]

    def get_nonpolar_score(self):
        return self.cluster.nonpolar_score[self.index]

    def get_total_score(self):
        return self.cluster.total_score[self.index]

    def get_contact(self):
        return self.contact

    def set_contact(self, contact):
        self.contact = contact



    target.get_pocket = get_pocket
    target.set_cluster = set_cluster
    target.get_cluster = get_cluster
    target.get_polar_score = get_polar_score
    target.get_nonpolar_score = get_nonpolar_score
    target.get_contact = get_contact
    target.set_contact = set_contact


def update_residue_method(target):
    def get_alphas(self):
        return self.atoms

    def get_polar_score(self):
        return np.sum([alpha.get_polar_score() for alpha in self.get_alphas()])

    target.get_polar_score = get_polar_score


    def get_nonpolar_score(self):
        return np.sum([alpha.get_nonpolar_score() for alpha in self.get_alphas()])

    def get_total_score(self):
        return np.sum([alpha.get_total_score() for alpha in self.get_alphas()])

    def set_cluster(self, cluster):
        self.cluster = cluster

    target.set_cluster = set_cluster

    def parent(self):
        return self.cluster

    def get_cluster(self):
        return self.cluster

    def get_contact(self):
        return np.any([alpha.get_contact() for alpha in self.get_alphas()])

    def get_alpha_index(self):
        return [atom.index for atom in self.get_alphas]

    def get_lining_atoms(self):
        return self.cluster._get_pocket_lining_atoms(self)
    target.get_lining_atoms = get_lining_atoms

    def get_lining_residues(self):
        return self.cluster._get_pocket_lining_residues(self)
    target.get_lining_residues = get_lining_residues

    target.get_nonpolar_score = get_nonpolar_score
    target.get_contact = get_contact
    target.get_alphas = get_alphas
    target.parent = parent
    target.get_total_score = get_total_score
