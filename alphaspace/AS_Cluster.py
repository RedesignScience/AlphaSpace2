"""
This file contains the container classes for cluster based objects.

AS_Data: A numpy.ndarray inheritance that stores and enumerates data
AS_Pocket: A mask container for pocket information, gets data from AS_Data
AS_AlphaAtom: A mast container for alpha atom, gets data from AS_Data
AS_BetaATOM: same as alpha atom
"""
import math

import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster

ASDATA_idx = 0
ASDATA_snapshot_idx = 1
ASDATA_x = 2
ASDATA_y = 3
ASDATA_z = 4
ASDATA_lining_atom_idx_1 = 5
ASDATA_lining_atom_idx_2 = 6
ASDATA_lining_atom_idx_3 = 7
ASDATA_lining_atom_idx_4 = 8
ASDATA_polar_space = 9
ASDATA_nonpolar_space = 10
ASDATA_is_active = 11
ASDATA_is_contact = 12
ASDATA_pocket_idx = 13
ASDATA_radii = 14
ASDATA_closest_binder_atom_idx = 15
ASDATA_closest_binder_atom_distance = 16
ASDATA_total_lining_atom_asa = 17


class AS_Data(np.ndarray):
    """
    Container object inherited from numpy array object, you can access information directly here.
    Most of AS_Data runs under the hood so you can safely ignore this part until needed.

    Column, Content,   default value
    0       idx
    1       snapshot_idx
    2       x
    3       y
    4       z
    5       lining_atom_idx_1
    6       lining_atom_idx_1
    7       lining_atom_idx_1
    8       lining_atom_idx_1
    9       polar_space        0
    10      nonpolar_space     0
    11      is_active          1
    12      is_contact         0
    13      pocket_idx         0
    14      radii
    15      closest binder atom idx   -1
    16      closest binder distance   0
    17      total lining atom asa
    """

    # noinspection PyArgumentList
    def __new__(cls, array_data, parent_structure=None):
        obj = super(AS_Data, cls).__new__(cls, buffer=array_data, shape=array_data.shape)
        obj.parent_structure = parent_structure
        return obj

    def idx(self, idx):
        """
        Get the absolute index of an item
        :param idx: int
        :return: int
        """
        return self[np.array(idx), 0]

    def snapshot_idx(self, idx: int):
        """
        Get the snapshot_index of an row
        :param idx: int
        :return: int
        """
        return self[np.array(idx), 1]

    @property
    def snapshots_idx(self):
        """
        Get all the snapshots in this data
        :return: np.ndarray
        """
        return np.sort(np.unique(self[:, 1])).astype(int)

    def snapshot(self, snapshot_idx):
        """
        Get all entries with in the given snapshot
        :param snapshot_idx: int
        :return: np.ndarray
        """
        return self[self.snapshot_alpha_idx(snapshot_idx)]

    def xyz(self, idx=None):
        """
        Get the coordinate of an entry
        :param idx: int or numpy.array or None
        :return: array of 3 or 3*N or all xyz
        """
        if idx is not None:
            return self[np.array(idx), 2:5]
        else:
            return self[:, 2:5]

    def get_polar_space(self, idx):
        """
        get the polar space value of a given entry
        :param idx:  int or numpy.array
        :return: polar_space
        :rtype: alphaspace.AS_Cluster.AS_Data
        """
        return self[np.array(idx), 9]

    def get_nonpolar_space(self, idx):

        return self[np.array(idx), 10]

    def lining_atoms_idx(self, idx):
        return np.array(self[np.array(idx, dtype=int), 5:9], dtype=int)

    def is_active(self, idx):
        return self[np.array(idx), 11]

    def is_contact(self, idx):
        return self[np.array(idx), 12]

    def pocket(self, idx):
        return self[np.array(idx), 13].astype(int)

    def snapshot_pockets_idx(self, snapshot_idx=0):
        return np.array(self[self[:, 1] == snapshot_idx][:, 13], dtype=int)

    def snapshot_alpha_idx(self, snapshot_idx=0):
        return np.array(self[self[:, 1] == snapshot_idx][:, 0], dtype=int)

    def radii(self, idx):
        return self[np.array(idx), 14]

    def activate(self, idx):
        self[idx, 11] = 1

    def deactivate(self, idx):
        self[idx, 11] = 0

    @property
    def all_idx(self):
        return self[:, 0].astype(int)

    def get_closest_atom_idx(self, idx):
        return self[idx, 15].astype(int)

    def get_closest_atom_dist(self, idx):
        return self[idx, 16].astype(float)

    def set_closest_atom(self, idx, values):
        assert len(idx) == len(values)
        self[idx, 15] = values

    def lining_atom_asa(self, idx):
        return self[idx, 17]


class AS_AlphaAtom:
    """
    This is the alpha atom container, which is the constituents of the pocket object.
    This object only contains reference of alpha atom id, and all info is stored in the array you can access by methods.
    """

    def __init__(self, idx, parent_structure, parent_pocket=None):
        self._idx = idx
        self.parent_structure = parent_structure
        self.parent_pocket = parent_pocket

    def __repr__(self):
        return "AlphaAtom #{}".format(self._idx)

    def __bool__(self):
        return self.is_active

    def __sub__(self, other):
        return np.linalg.norm(self.xyz - other.xyz)

    @property
    def _data(self) -> AS_Data:
        return self.parent_structure._data

    @property
    def is_active(self):
        return bool(self._data.is_active(self._idx))

    @property
    def idx(self):
        return self._idx

    @property
    def snapshot_idx(self):
        return self._data.snapshot(self._idx)

    @property
    def xyz(self):
        return self._data.xyz(self._idx)

    @property
    def polar_space(self):
        return self._data.get_polar_space(self._idx)

    @property
    def nonpolar_space(self):
        return self._data.get_nonpolar_space(self._idx)

    @property
    def space(self):
        return self.polar_space + self.nonpolar_space

    @property
    def lining_atoms_idx(self):
        return self._data.lining_atoms_idx(self._idx)

    @property
    def is_contact(self) -> bool:
        return bool(self._data.is_contact(self._idx))

    @property
    def lining_atom_asa(self):
        return float(self._data.lining_atoms_idx(self._idx))

    @property
    def is_solvated(self):
        return self.lining_atom_asa > 0

    @property
    def closest_atom_idx(self):
        return self._data.get_closest_atom_idx(self._idx)

    @property
    def closest_atom_dist(self):
        return self._data.get_closest_atom_dist(self._idx)


# noinspection PyTypeChecker
class AS_Pocket:
    """
    This is the pocket container for the topological information of a pocket
    A pocket is a collection of alpha atoms clustered by average linkage
    """

    def __init__(self, alpha_idx, snapshot_idx, idx, parent_structure):
        self._idx = idx
        self._snapshot_idx = snapshot_idx
        self._alpha_idx = alpha_idx
        self.parent_structure = parent_structure

        self._contact_pockets = set()
        self._lining_atoms_idx = None
        self._reordered_index = None

        self._connected = False
        self._is_core_d = False

    @property
    def snapshot_idx(self):
        return self._snapshot_idx

    def set_contact(self, other):
        """
        Mark two pockets as contact and link them together.
        :param other: the other pocket

        """
        self._contact_pockets.add(other.index)
        other.contact_pockets.add(self.index)

    def __ne__(self, other):
        return not self.__eq__(other)

    def __repr__(self):
        return "<Pocket {} in snapshot {} with space {} and {}% occupied>".format(self.reordered_index,
                                                                                  self._snapshot_idx,
                                                                                  math.ceil(self.space),
                                                                                  math.ceil(self.occupancy * 100))

    @property
    def index(self):
        return int(self._idx)

    @property
    def universe(self):
        return self.parent_structure.parent.universe

    @property
    def reordered_index(self):
        """
        This returns the index of the pocket in the snapshot AFTER ranking them by space, and is mutable.
        If the reordered_index hasn't been set, original index is returned.
        :return: int
        """
        if self._reordered_index is None:
            return self.index
        else:
            return int(self._reordered_index)

    def __hash__(self):
        return hash((self._snapshot_idx, self._idx))

    def __eq__(self, other):
        if (self._idx, self._snapshot_idx) == (other._idx, other._snapshot_idx):
            return True
        else:
            return False

    def __len__(self) -> int:
        return len(self._alpha_idx)

    def activate(self):
        """
        activate this pocket, and all the member alpha atoms in the pocket.
        """
        self._data.activate(self.alpha_idx)

    def deactivate(self):
        """
        deactivate this pocket, and all the member alpha atoms in the pocket.
        """
        self._data.deactivate(self.alpha_idx)

    @property
    def config(self):
        return self.parent_structure.config

    @property
    def alphas(self):
        """
        Generator of all alpha atoms in the pocket
        :return: alp
        """
        for i in self.alpha_idx:
            yield AS_AlphaAtom(i, self.parent_structure, self)

    @property
    def alpha_idx(self):
        return np.array(self._alpha_idx, dtype=int)

    @property
    def _data(self):
        return self.parent_structure._data

    @property
    def is_active(self, anyone=True):
        if anyone:
            return self._data.is_active(self.alpha_idx).any()
        else:
            return self._data.is_active(self.alpha_idx).all()

    @property
    def xyz(self) -> np.ndarray:
        return self._data.xyz(self.alpha_idx)

    @property
    def lining_atoms_idx(self):
        if self._lining_atoms_idx is None:
            self._lining_atoms_idx = np.unique(self._data.lining_atoms_idx(self.alpha_idx))
        return self._lining_atoms_idx

    @property
    def lining_atoms_idx_redu(self):
        return self._data.lining_atoms_idx(self.alpha_idx)

    @property
    def lining_atoms(self):
        for idx in self.lining_atoms_idx:
            yield self.parent_structure.top.atom(idx)

    @property
    def lining_atoms_xyz(self):
        return self.parent_structure.traj.xyz[self._snapshot_idx, self.lining_atoms_idx, :]

    @property
    def lining_atoms_centroid(self):
        return np.average(self.lining_atoms_xyz, axis=0)

    @property
    def lining_residues_idx(self):
        return np.unique([atom.residue.index for atom in self.lining_atoms])

    @property
    def lining_residues(self):
        residue_idx = self.lining_residues_idx
        for idx in residue_idx:
            yield self.parent_structure.top.residue(idx)

    @property
    def color(self):
        return self.parent_structure.config.color_name(self._idx)

    @property
    def is_contact(self):
        return np.any(self.contacts)

    @property
    def occupancy(self):
        if self.is_contact:
            return float(self.get_space(contact_only=True, space_type='ALL') / self.get_space(contact_only=False,
                                                                                              space_type='ALL'))
        else:
            return 0.0

    @property
    def contacts(self) -> np.ndarray:
        return np.array(self._data.is_active(self.alpha_idx), dtype=int)

    def get_spaces(self, space_type='ALL', contact_only=False):
        if space_type.upper() == 'POLAR':
            spaces = self.get_polar_spaces()
        elif space_type.upper() == 'NONPOLAR':
            spaces = self.get_nonpolar_spaces()
        elif space_type.upper() == 'ALL':
            spaces = self.get_total_spaces()
        else:
            raise Exception("Use ALL, POLAR, NONPOLAR as space_type")
        if contact_only:
            return self.contacts * spaces
        else:
            return spaces

    def get_space(self, space_type='ALL', contact_only=False):
        return float(np.sum(self.get_spaces(space_type, contact_only)))

    def get_polar_space(self) -> float:
        return float(np.sum(self._data.get_polar_space(self.alpha_idx)))

    def get_polar_spaces(self) -> np.ndarray:
        return self._data.get_polar_space(self.alpha_idx)

    def get_nonpolar_space(self) -> float:
        return float(np.sum(self._data.get_nonpolar_space(self.alpha_idx)))

    def get_nonpolar_spaces(self) -> np.ndarray:
        return self._data.get_nonpolar_space(self.alpha_idx)

    def get_total_space(self):
        return self.get_polar_space() + self.get_nonpolar_space()

    def get_total_spaces(self):
        return self.get_polar_spaces() + self.get_nonpolar_spaces()

    @property
    def polar_space(self):
        return self.get_polar_space()

    @property
    def nonpolar_space(self):
        return self.get_nonpolar_space()

    @property
    def space(self):
        return self.get_space()

    @property
    def lining_atom_asa(self) -> np.ndarray:
        return self._data.lining_atom_asa(self.lining_atoms_idx)

    @property
    def betas(self):
        if len(self) > 1:
            zmat = linkage(self.xyz, method='average')
            beta_cluster_idx = fcluster(zmat, self.config.beta_clust_cutoff / 10, criterion='distance') - 1
            beta_list = [[] for _ in range(max(beta_cluster_idx) + 1)]
            for i, c_i in enumerate(beta_cluster_idx):
                beta_list[c_i].append(i)
            for beta in beta_list:
                yield AS_BetaAtom(alpha_idx_in_pocket=beta, pocket=self)
        else:
            yield AS_BetaAtom([0], self)

    @property
    def centroid(self):
        return np.average(self.xyz, axis=0)

    @property
    def core_aux_minor(self):
        if self.space < self.config.aux_cutoff:
            return 'minor'
        elif self.space < self.config.core_cutoff:
            return 'aux'
        else:
            return 'core'

    def union(self, other) -> set:
        return set(self.lining_atoms_idx).union(set(other.lining_atoms_idx))

    def intersection(self, other) -> set:
        return set(self.lining_atoms_idx).intersection(set(other.lining_atoms_idx))

    def jaccard_similarity(self, other) -> float:
        return float(len(self.intersection(other))) / len(self.union(other))


class AS_BetaAtom:
    def __init__(self, alpha_idx_in_pocket: list, pocket: AS_Pocket):
        self._pocket = pocket
        self._alpha_idx_in_pocket = alpha_idx_in_pocket
        self.alpha_idx = pocket.alpha_idx[alpha_idx_in_pocket]

    @property
    def xyz(self):
        return self.pocket.xyz[self._alpha_idx_in_pocket]

    @property
    def centroid(self):
        return np.average(self.xyz, axis=0)

    @property
    def alphas(self):
        for alpha_idx in self.alpha_idx:
            yield AS_AlphaAtom(idx=alpha_idx, parent_structure=self.pocket.parent_structure, parent_pocket=self.pocket)

    @property
    def pocket(self):
        return self._pocket

    @property
    def space(self):
        return np.sum([alpha.space for alpha in self.alphas])


class AS_D_Pocket:
    def __init__(self, parent_structure, pocket_indices=None):
        """
        Initialize a mask for information storage of a dpocket, this is user accessible and can only be read from public
         methods
        :param parent_structure: AS_Universe, parent universe
        """
        self.structure = parent_structure
        if pocket_indices is not None:
            self._pocket_indices = pocket_indices  # list of tuple, (snapshot_idx, pocket_idx)
        self.index = -1  # index of d-pocket in all d-pocket
        self._active = True
        self._data = self.structure._data

    @property
    def is_active(self):
        return self._active

    def activate(self):
        self._active = True

    def deactivate(self):
        self._active = False

    def __bool__(self):
        return self._active

    @property
    def _pocket_set(self):
        return set(self._pocket_indices)

    def __getitem__(self, key):
        """
        :param key:
        :return:
        """

    def _index_to_pocket(self, index: tuple):
        return self.structure.cluster(snapshot_idx=index[0]).pocket(index[1])

    def _pocket_to_index(self, pocket) -> tuple:
        return pocket.cluster.snapshot_index, pocket.index

    def __contains__(self, item) -> bool:
        if type(item) is not tuple:
            item = self._pocket_to_index(item)
        return any(p == item for p in self._pocket_indices)

    @property
    def __iter__(self) -> object:
        for idx in self._pocket_indices:
            yield self._index_to_pocket(idx)

    @property
    def pocket_indices(self) -> tuple:
        for i in self._pocket_indices:
            yield i

    @property
    def pockets(self):
        for pocket_idx in self.pocket_indices:
            yield self._index_to_pocket(pocket_idx)

    @property
    def pockets_by_snapshots(self):
        snapshot_dict = {}
        for ss_idx, pkt_idx in self._pocket_indices:
            if ss_idx not in snapshot_dict:
                snapshot_dict[ss_idx] = []
            snapshot_dict[ss_idx].append(self._index_to_pocket((ss_idx, pkt_idx)))
        return snapshot_dict

    def add(self, item):
        """
        add a pocket to this dpocket
        :param item: object, AS_Pocket
        :return: bool, if the pocket isn't already in the dpocket list
        """
        if type(item) is not tuple:
            item = self._pocket_to_index(item)
        if item in self:
            raise Exception('Pocket {} already exist in Frame {}'.format(item[1], item[0]))

    def pop(self, key):
        """
        pop a pocket
        :param key: int
        :return: AS_Pocket
        """
        return self._index_to_pocket(self._pocket_indices.pop(key))

    @property
    def n_pockets(self) -> int:
        """
        total number of pockets
        :return: int
        """
        return len(self._pocket_indices)

    @property
    def n_snapshots(self) -> int:
        """
        total number of covered snapshot
        :return:
        """
        return len(self.snapshots)

    @property
    def snapshots(self) -> list:
        """
        Return snapshot index for all the pockets in this d-pocket
        :return: list of int
        """
        return sorted(list(set([i[1] for i in self._pocket_indices])))

    @property
    def pocket_spaces(self) -> tuple:
        """
        :return: iter of tuple (polar, nonpolar)
        """
        for pocket in self.pockets:
            yield (pocket.get_polar_space(), pocket.get_nonpolar_space())

    @property
    def pocket_xyz(self) -> np.array:
        for pocket in self.pockets:
            yield pocket.get_centroid()
