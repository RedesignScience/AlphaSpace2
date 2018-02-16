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

import alphaspace

class AS_Snapshot(np.ndarray):
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
    def __new__(cls, array_data):
        obj = super(AS_Snapshot, cls).__new__(cls, buffer=array_data, shape=array_data.shape)
        return obj

    def idx(self, idx):
        """
        Get the absolute index of an item
        :param idx: int
        :return: int
        """
        return self[np.array(idx), 0]

    def snapshot_idx(self):
        """
        Get the snapshot_index
        :param idx: int
        :return: int
        """
        return self[0, 1]

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

    def is_active(self, idx=None):
        if idx is None:
            return self[:, 11]
        else:
            return self[np.array(idx), 11]

    def is_contact(self, idx):
        return self[np.array(idx), 12]

    def pocket(self, idx):
        return self[np.array(idx), 13].astype(int)

    def snapshot_pockets_idx(self):
        return np.array(self[:, 13], dtype=int)

    def snapshot_alpha_idx(self):
        return np.array(self[:, 0], dtype=int)

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


class AS_Data():
    def __init__(self, data_dict, parent_structure=None):
        self.parent_structure = parent_structure
        self.data_dict = data_dict

    def update(self, data_dict):
        self.data_dict.update(data_dict)

    def __getitem__(self, i):
        try:
            return self.data_dict[i]
        except IndexError:
            print("No snapshot {}".format(i))
            return False

    @property
    def snapshots_idx(self):
        """
        Get all the snapshots in this data
        :return: np.ndarray
        """
        return sorted(self.data_dict.keys)

    def snapshot(self, snapshot_idx):
        """
        Get all entries with in the given snapshot
        :param snapshot_idx: int
        :return: np.ndarray
        """
        return self[snapshot_idx]

    def xyz(self, alpha_idx, snapshot_idx):
        """
        Get the coordinate of an entry
        :param alpha_idx: int or numpy.array or None
        :return: array of 3 or 3*N or all xyz
        """

        return self[snapshot_idx].xyz(alpha_idx)

    def get_polar_space(self, alpha_idx, snapshot_idx):
        """
        get the polar space value of a given entry
        :param idx:  int or numpy.array
        :return: polar_space
        :rtype: float
        """
        return self[snapshot_idx].get_polar_space(alpha_idx)

    def get_nonpolar_space(self, alpha_idx, snapshot_idx):
        return self[snapshot_idx].get_nonpolar_space(alpha_idx)

    def lining_atoms_idx(self, alpha_idx, snapshot_idx):
        return self[snapshot_idx].lining_atoms_idx(alpha_idx)

    def is_active(self, alpha_idx, snapshot_idx):
        return self[snapshot_idx].is_active(alpha_idx)

    def is_contact(self, alpha_idx, snapshot_idx):
        return self[snapshot_idx].xyz(alpha_idx)

    def pocket(self, alpha_idx, snapshot_idx):
        return self[snapshot_idx].pocket(alpha_idx)

    def snapshot_alpha_idx(self, snapshot_idx=0):
        return self[snapshot_idx].snapshot_alpha_idx()

    def radii(self, alpha_idx, snapshot_idx):
        return self[snapshot_idx].radii(alpha_idx)

    def activate(self, alpha_idx, snapshot_idx):
        self[snapshot_idx].activate(alpha_idx)

    def deactivate(self, alpha_idx, snapshot_idx):
        self[snapshot_idx].deactivate(alpha_idx)

    def get_closest_atom_idx(self, alpha_idx, snapshot_idx):
        return self[snapshot_idx].get_closest_atom_idx(alpha_idx)

    def get_closest_atom_dist(self, alpha_idx, snapshot_idx):
        return self[snapshot_idx].get_closest_atom_idx(alpha_idx)

    def set_closest_atom(self, alpha_idx, values, snapshot_idx):
        assert len(alpha_idx) == len(values)
        self[snapshot_idx].set_closest_atom(alpha_idx, values)

    def lining_atom_asa(self, alpha_idx, snapshot_idx):
        return self[snapshot_idx].lining_atom_asa(alpha_idx)


class AS_AlphaAtom:
    """
    This is the alpha atom container, which is the constituents of the pocket object.

    This object only contains reference of alpha atom id, and all info is stored in the array you can access by methods.
    """

    _ngl_radius = 0.5

    def __init__(self, idx, snapshot_idx, parent_structure, parent_pocket=None):
        """

        Parameters
        ----------
        idx: int
        snapshot_idx: int
        parent_structure : AS_Struct
        parent_pocket: AS_Pocket or None
        """
        self._idx = idx
        self._snapshot_idx = snapshot_idx
        self.parent_structure = parent_structure
        self.parent_pocket = parent_pocket

        self._ngl_component_idx = None

    def __repr__(self):
        return "AlphaAtom #{} in Snapshot {}".format(self._idx, self.snapshot_idx)

    def __bool__(self):
        return self.is_active

    def __sub__(self, other):
        """
        Use subtraction to get the distance between two alpha atoms

        Parameters
        ----------
        other : AS_AlphaAtom

        Returns
        -------
        distance : float
        """
        return np.linalg.norm(self.xyz - other.xyz)

    @property
    def _data(self) -> AS_Snapshot:
        return self.parent_structure._data[self.snapshot_idx]

    @property
    def is_active(self):
        """
        See if this alpha atom is activated, default is based on screening threshold.

        Returns
        -------
        bool

        """
        return bool(self._data.is_active(self.idx))

    @property
    def idx(self):
        """
        Get the index

        Returns
        -------
        index : int
        """
        return self._idx

    @property
    def snapshot_idx(self):
        """
        Get the parent snapshot index

        Returns
        -------

        snapshot_idx: float

        """
        return self._snapshot_idx

    @property
    def xyz(self):
        """

        Returns
        -------

        xyz : np.ndarray
            size 3

        """
        return self._data.xyz(self._idx)

    @property
    def centroid(self):
        return self.xyz

    @property
    def polar_space(self):
        """

        Returns
        -------
        float

        """
        return self._data.get_polar_space(self._idx)

    @property
    def nonpolar_space(self):
        """

        Returns
        -------
        float

        """
        return self._data.get_nonpolar_space(self._idx)

    @property
    def space(self):
        """

        Returns
        -------
        space : float


        """
        return self.polar_space + self.nonpolar_space

    @property
    def lining_atoms_idx(self):
        """
        Get the atom index of all the lining atoms. Each alpha atom has four lining atoms.


        Returns
        -------

        indices : np.ndarray
            size (4)

        """
        return self._data.lining_atoms_idx(self._idx)

    @property
    def is_contact(self) -> bool:
        """
        Check if it's in contact with any binder atoms

        Returns
        -------

        bool

        """
        return bool(self._data.is_contact(self._idx))

    @property
    def lining_atom_asa(self):
        """
        Get total lining atom solvent accessible area

        Returns
        -------
        ASA : float
        """

        return float(self._data.lining_atoms_idx(self._idx))

    @property
    def is_solvated(self):
        """
        If the lining atom can be solvated,
        True if lining_atom_asa > 0

        Returns
        -------
        bool

        """
        return self.lining_atom_asa > 0

    @property
    def closest_atom_idx(self):
        """

        Returns
        -------
        index : np.array

        """
        return self._data.get_closest_atom_idx(self._idx)

    @property
    def closest_atom_dist(self):
        """

        Returns
        -------
        distance to closest atom : float
        """
        return self._data.get_closest_atom_dist(self._idx)


# noinspection PyTypeChecker
class AS_Pocket:
    """
    This is the pocket container for the topological information of a pocket
    A pocket is a collection of alpha atoms clustered by average linkage


    Initialize pocket object by calling .pockets or .pocket method in AS_Universe.

    Examples
    --------

    If you want to see pockets from snapshot # 3

    >>>   import alphaspace
    >>>   universe = alphaspace.AS_Universe()
    >>>   for pocket in universe.pockets(snapshot_idx = 3, active_only = True):
    >>>      print(pocket)

    The pocket object is a mask container linked to AS_Data,
    you can use the build in methods to access the data such as contact, occupancy, space or geometric info.


    """

    _ngl_radius = 2.0

    def __init__(self, alpha_idx, snapshot_idx, pocket_idx, parent_structure):
        self._idx = pocket_idx
        self._snapshot_idx = snapshot_idx
        self._alpha_idx = alpha_idx
        self.parent_structure = parent_structure
        self._contact_pockets = set()
        self._lining_atoms_idx = None
        self._reordered_index = None
        self._connected = False
        self._is_core_d = False
        self._betas = None

        self._ngl_component_idx = None


    @property
    def snapshot_idx(self):
        """

        Returns
        -------

        snapshot_idx : int
            The index of snapshot this pocket in located

        """
        return self._snapshot_idx

    def set_contact(self, other):
        """
        Mark two pockets as contact and link them together.

        Parameters
        ----------

        other : AS_Pocket
            The other pocket to be linked

        """
        self._contact_pockets.add(other.index)
        other._contact_pockets.add(self.index)

    @property
    def index(self):
        """
        Gives the absolute index of this pocket in the snapshot,
        this index is used fo pocket look up.

        Returns
        -------

        int

        """
        return int(self._idx)

    @property
    def universe(self):
        """
        Traceback to the parent universe

        Returns
        -------

        universe : AS_Universe
            parent universe

        """
        return self.parent_structure.universe

    @property
    def reordered_index(self):
        """
        This returns the index of the pocket in the snapshot AFTER ranking them by space, and is mutable.
        If the reordered_index hasn't been set, original index is returned.

        Returns
        -------

        int

        """
        if self._reordered_index is None:
            return self.index
        else:
            return int(self._reordered_index)

    def __repr__(self):
        return "<Pocket {} in snapshot {} with space {} and {}% occupied>".format(self.reordered_index,
                                                                                  self._snapshot_idx,
                                                                                  math.ceil(self.space),
                                                                                  math.ceil(self.occupancy * 100))

    def __hash__(self):
        return hash((self._snapshot_idx, self._idx))

    def __eq__(self, other):
        if (self._idx, self._snapshot_idx) == (other._idx, other._snapshot_idx):
            return True
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def __len__(self) -> int:
        """
        Number of alpha atoms

        Returns
        -------
        int
        """
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
        """
        Get the configuration container

        Returns
        -------

        AS_Config

        """
        return self.parent_structure.config

    @property
    def alphas(self):
        """
        Generator of all alpha atoms in the pocket

        Returns
        -------

        iterable : AS_AlphaAtom

        """
        for i in self.alpha_idx:
            yield AS_AlphaAtom(i,snapshot_idx=self.snapshot_idx,parent_structure=self.parent_structure,parent_pocket=self)

    @property
    def alpha_idx(self):
        """
        Get a numpy array of child alpha atom indices

        Returns
        -------

        indicies : np.ndarray

        """
        return np.array(self._alpha_idx, dtype=int)

    @property
    def _data(self):
        return self.parent_structure._data[self.snapshot_idx]

    @property
    def is_active(self, anyone=True):
        """
        Toggle for if pocket is active

        Parameters
        ----------

        anyone : bool
            default if True, any active alpha atom in the pocket will toggle the pocket as active.
            false requires all alphas to be active.

        Returns
        -------

        bool

        """
        if anyone:
            return self._data.is_active(self.alpha_idx).any()
        else:
            return self._data.is_active(self.alpha_idx).all()

    @property
    def xyz(self) -> np.ndarray:
        """
        Get the coordinates for all alphas
        Returns
        -------

        xyz : np.ndarray
            shape (n, 3) for n alpha atoms

        See Also
        --------
        centroid : get the centroid of all alpha atoms in the pocket.
        """
        return self._data.xyz(self.alpha_idx)

    @property
    def lining_atoms_idx(self):
        if self._lining_atoms_idx is None:
            self._lining_atoms_idx = np.unique(self._data.lining_atoms_idx(self.alpha_idx))
        return self._lining_atoms_idx

    @property
    def lining_atoms_vector(self):

        vector = np.zeros(self.parent_structure.n_atoms)

        np.put(vector, self.lining_atoms_idx, 1)
        return vector

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
        """
        Get the color for visualization

        Returns
        -------

        color_name : str
        """
        return self.parent_structure.config.color_name(self._idx)

    @property
    def is_contact(self):
        """

        Returns
        -------

        bool
        """
        return np.any(self.contacts)

    @property
    def occupancy(self):
        """
        Get the percentage of alpha atoms in contact with binder.

        Returns
        -------
        occupancy : float
            0 if there is no binder

        """
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
        # noinspection PyUnresolvedReferences
        """
        Iterate over beta atoms in this pocket.

        Yields
        ------
        beta_atom : AS_BetaAtom

        See Also:
        ---------
        AS_BetaAtom

        """
        if self._betas is None:
            if len(self) > 1:
                zmat = linkage(self.xyz, method='average')
                beta_cluster_idx = fcluster(zmat, 1.6 / 10, criterion='distance') - 1
                beta_list = [[] for _ in range(max(beta_cluster_idx) + 1)]
                for i, c_i in enumerate(beta_cluster_idx):
                    beta_list[c_i].append(i)
                self._betas = [AS_BetaAtom(alpha_idx_in_pocket=beta, pocket=self) for beta in beta_list]
            else:
                self._betas = [AS_BetaAtom([0], self)]
        return iter(self._betas)


    @property
    def score(self):
        """
        Calculate the ligandibility for this pocket, it's the sum of all beta atom scores

        Returns
        -------
        score: float

        """

        return np.sum([beta.score for beta in self.betas])

    @property
    def centroid(self):
        """
        Calculate the centroid of all alpha atoms in this pocket.

        Returns
        -------
        centroid : np.ndarray
            shape : (3)
        """

        return np.mean(self.xyz, axis=0)

    @property
    def core_aux_minor(self):
        if self.space < self.config.aux_cutoff:
            return 'minor'
        elif self.space < self.config.core_cutoff:
            return 'aux'
        else:
            return 'core'

    def union(self, other) -> set:
        """

        Parameters
        ----------
        other

        Returns
        -------

        """
        return set(self.lining_atoms_idx).union(set(other.lining_atoms_idx))

    def intersection(self, other) -> set:
        return set(self.lining_atoms_idx).intersection(set(other.lining_atoms_idx))

    def jaccard_similarity(self, other) -> float:
        return float(len(self.intersection(other))) / len(self.union(other))


class AS_BetaAtom:
    """
    This is the container for Beta atom, which is simply a collection of alpha atoms.

    It belongs to the AS_Pocket object.
    """

    _ngl_radius = 1.0


    def __init__(self, alpha_idx_in_pocket, pocket: AS_Pocket):
        """
        AS_BetaAtom are automatically generated from 'AS_Pocket.betas' iterator.

        Parameters
        ----------
        alpha_idx_in_pocket : list
        pocket : AS_Pocket
        """
        self._pocket = pocket
        self._alpha_idx_in_pocket = alpha_idx_in_pocket
        self.alpha_idx = pocket.alpha_idx[alpha_idx_in_pocket]

        self.prb_element = []
        self._vina_score = None

        self._ngl_component_idx = None

    def __repr__(self):
        return "Beta atom with {} space".format(self.space)

    @property
    def xyz(self):
        """
        Gets a list of all alpha atom coordinates

        Returns
        -------

        coordinates : np.ndarray
            shape = (n,3) for n alpha atoms

        """

        return self.pocket.xyz[self._alpha_idx_in_pocket]

    @property
    def centroid(self):
        """
        Gets the centroid of this beta atom


        Returns
        -------

        centroid coordinate : np.ndarray
            shape = (3,)
        """
        return np.mean(self.xyz, axis=0)

    @property
    def alphas(self):
        """
        Iterate over member alpha atoms

        Yields
        ------
        AlphaAtom : AS_AlphaAtom

        """
        for alpha_idx in self.alpha_idx:
            yield AS_AlphaAtom(idx=alpha_idx, snapshot_idx=self._pocket.snapshot_idx ,parent_structure=self.pocket.parent_structure, parent_pocket=self.pocket)

    @property
    def pocket(self):
        """
        Get the parent pocket

        Returns
        -------
        parent_pocket : AS_Pocket

        """

        return self._pocket

    @property
    def space(self):
        """
        Total space of all alphas

        Returns
        -------

        total_space : float

        """
        return np.sum([alpha.space for alpha in self.alphas])

    @property
    def vina_scores(self):
        """
        Return all terms of the Vina Scores.

        Returns
        -------
        vina_scores : np.ndarray
            shape : (9 , 6)
            9 probe element by 6 scores
            PROBE_TYPE = ['C', 'Br', 'F', 'Cl', 'I', 'OA', 'SA', 'N', 'P']
            scores ={'total', "gauss_1","gauss_2", "repulsion", "hydrophobic", "Hydrogen"}


        """


        if self._vina_score is None:
            raise Exception('No Vina Score Calculated')
        else:
            return self._vina_score

    @property
    def best_probe_type(self):

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

        return alphaspace.best_probe_type(self)



    @property
    def score(self):
        """
        Get the score of this beta atom, which is the lowest vina score in all 9 probes.

        Returns
        -------
        float

        """
        return np.min(self.vina_scores[:,0])


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

        self._ngl_component_idx = None

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


class AS_community:
    def __init__(self, core = None, aux = None):
        self._core = core if core is not None else []
        self._aux = aux if aux is not None else []
        self._lining_atoms = None

    def connected_core(self, pocket):
        from .AS_Funct import is_pocket_connected
        for c in self._core:
            if is_pocket_connected(c, pocket):
                return True
        return False

    def __len__(self):
        return len(self.pockets)

    @property
    def pockets(self):
        return self._core + self._aux

    @property
    def core(self):
        return self._core
    @property
    def aux(self):
        return self._aux

    @property
    def lining_atoms(self):
        if self._lining_atoms is None:
            lining_atoms = np.hstack([pocket.lining_atoms_idx for pocket in self.pockets])
            self._lining_atoms = np.unique(lining_atoms)
        return self._lining_atoms


    def similarity(self,community):

        return len(np.intersect1d(self.lining_atoms,community.lining_atoms,assume_unique=True))/len(
            np.unique(np.hstack([self.lining_atoms, community.lining_atoms])))

    def similarity_atoms(self,lining_atoms):
        return len(np.intersect1d(self.lining_atoms, np.array(lining_atoms), assume_unique=True)) / len(
            np.unique(np.hstack([self.lining_atoms, np.array(lining_atoms)])))

