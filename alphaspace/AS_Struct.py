"""
AS_Struct is the container for structures in the universe

AS_Universe:
            : receprot_struct
            : binder_struct
            : pockets and other virtual elements
"""

from collections import defaultdict

import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial import Voronoi, Delaunay
from scipy.spatial.distance import squareform
from .AS_Cluster import AS_D_Pocket, AS_Data, AS_Pocket
from .AS_Funct import getTetrahedronVolume, getSASA, getIfContact


class AS_Structure:
    def __init__(self, trajectory, structure_type: int = 2, parent=None):

        """
        Container for structure trajectory and topology in a AS_Session
        :param trajectory: MDtraj trajectory object
        :param structure_type: int 0 for receptor, 1 for binder, 2 for other
        """

        # 0 for receptor, 1 for binder, 2 for unassigned
        self.structure_type = structure_type
        self.trajectory = trajectory
        self.parent = parent
        self.universe = parent
        # self.contact_cluster = [[None for i in range(self.n_residues)] for j in range(self.n_frames)]
        self._data = None
        self._pockets_alpha_idx = {}
        self._pockets = {}

    def __bool__(self):
        return True

    def __repr__(self):
        return "{} Structure with {} frames, {} residues, {} atoms".format(
            ['Receptor', 'Binder', 'Misc.'][self.structure_type], self.n_frames, self.n_residues, self.n_atoms)

    def __len__(self):
        """
        Returns number of frames
        :return: int
        """
        return self.n_frames

    @property
    def config(self):
        return self.parent.config

    @property
    def traj(self):
        """
        :return: trajectory of this structure
        """
        return self.trajectory

    @property
    def n_atoms(self):
        """
        Get total number of atoms
        :return: int
        """
        return self.trajectory.n_atoms

    @property
    def n_frames(self):
        """
        Get the total number of frames.
        Same as n_snapshots
        :return: int
        """
        return self.trajectory.n_frames

    @property
    def n_snapshots(self):
        """
        Get the total number of snapshots.
        Same as n_frames
        :return: int
        """
        return self.trajectory.n_frames

    @property
    def n_residues(self):
        """
        Get the total number of residues in the structure
        :return:
        """
        return self.topology.n_residues

    @property
    def topology(self):
        """
        :return: topology of this structure
        """
        return self.trajectory.topology

    @property
    def top(self):
        """
        :return: topology of this structure
        """
        return self.trajectory.topology

    @property
    def residues(self):
        """
        Residue iterator
        :return: iter residue topology
        """
        return self.topology.residues

    @property
    def atoms(self):
        """
        Atom iterator
        :return: iter
        """
        return self.top._atoms
        #
        # for atom in self.top._atoms:
        #     yield atom

    @property
    def is_polar(self):
        """
        Returns an array of if the atom in the topology is a polar atom.
        :return: np.ndarray N of n atoms in the structure
        """
        return np.array([(str(atom.element) in ['nitrogen', 'oxygen', 'sulfur']) for atom in self.topology._atoms])

    def residue(self, idx):
        """
        Gives a residue with idx
        :param idx: int
        :return: Residue
        """
        return self.topology.residue(idx)

    def atom(self, idx):
        """
        Gives an atom with idx
        :param idx: int
        :return: object atom
        """
        return self.topology.atom(idx)

    def n_alphas(self, snapshot_idx=0, active_only=False):
        """
        return the number of alpha atoms in a particular snapshot
        :param snapshot_idx: int
        :param active_only: bool
        :return: int
        """

        if active_only:
            return len(np.where(self._data[snapshot_idx].is_active()))
        else:
            return len(self._data[snapshot_idx])

    def _gen_pockets(self):

        self._pockets_alpha_idx ={}

        for i in range(self.n_frames):
            pocket_snapshot_dict = self._data[i][:, [0, 13]]

            reversed_dict = defaultdict(list)
            for idx, p_idx in pocket_snapshot_dict:
                reversed_dict[p_idx].append(idx)
            self._pockets_alpha_idx[i] = reversed_dict

    def pockets(self, snapshot_idx=0):
        """
        Generate an iterator of the pockets in the given snapshot
        Parameters
        ----------
        snapshot_idx : int


        Returns
        -------
        iterable : alphaspace.AS_Pocket

        """
        if len(self._pockets_alpha_idx) == 0:
            self._gen_pockets()
        for pocket_idx, pocket_content in self._pockets_alpha_idx[snapshot_idx].items():
            yield AS_Pocket(pocket_content, snapshot_idx, pocket_idx, self)

    def pocket(self, pocket_idx, snapshot_idx=0):
        """
        Get a pocket by index and snapshot index
        None if it does not exist
        :param pocket_idx: int
        :param snapshot_idx: int
        :return: None or AS_Pocket
        """
        if self._pockets_alpha_idx:
            self._gen_pockets()
        if snapshot_idx in self._pockets_alpha_idx:
            if pocket_idx in self._pockets_alpha_idx[snapshot_idx]:
                return AS_Pocket(self._pockets_alpha_idx[snapshot_idx][pocket_idx], snapshot_idx, pocket_idx, self)
            else:
                return None
        else:
            return None

    def calculate_contact(self, snapshot_idx=None):
        """
        Calculate the contact index of the alpha cluster against the designated binder.
        The contact distance cutoff can be set in config
        :param snapshot_idx: int
        """

        if snapshot_idx is not None:

            snapshot_cluster_coord_matrix = self._data[snapshot_idx].xyz()
            binder_coords = self.universe.binder.trajectory.xyz[snapshot_idx]
            contact_alpha = getIfContact(snapshot_cluster_coord_matrix, binder_coords, self.config.hit_dist)[0]
            self._data[snapshot_idx][contact_alpha, 12] = 1

        else:
            for i in range(self.n_snapshots):
                self.calculate_contact(i)