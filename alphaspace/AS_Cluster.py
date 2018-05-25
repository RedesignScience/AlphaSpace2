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


class Alpha:
    """
    This is the alpha atom container, which is the constituents of the pocket object.

    This object only contains reference of alpha atom id, and all info is stored in the array you can access by methods.
    """

    def __init__(self, snapshot, index):
        """

        Parameters
        ----------
        snapshot: alpsphace.AS_Snapshot
        index: int
        """

        self.snapshot = snapshot
        self.index = index
        self.index = index

    def __repr__(self):
        return "AlphaAtom {}".format(self.index)

    def __sub__(self, other):
        """
        Use subtraction to get the distance between two alpha atoms

        Parameters
        ----------
        other : Alpha

        Returns
        -------
        distance : float
        """
        return np.linalg.norm(self.xyz - other.xyz)

    @property
    def xyz(self):
        """

        Returns
        -------

        xyz : np.ndarray
            size 3

        """
        return self.snapshot.alpha_xyz[self.index]

    @property
    def centroid(self):
        return self.xyz

    @property
    def space(self):
        """

        Returns
        -------
        space : float


        """
        return self.snapshot.alpha_space[self.index] * 1000

    @property
    def lining_atoms_idx(self):
        """
        Get the atom index of all the lining atoms. Each alpha atom has four lining atoms.


        Returns
        -------

        indices : np.ndarray
            size (4)

        """
        return self.snapshot.alpha_lining(self.index)

    @property
    def is_contact(self) -> bool:
        """
        Check if it's in contact with any binder atoms

        Returns
        -------

        bool

        """
        return bool(self.snapshot.alpha_contact[self.index])


class Beta:
    """
    This is the container for Beta atom, which is simply a collection of alpha atoms.

    It belongs to the AS_Pocket object.
    """

    def __init__(self, snapshot, index):
        self.snapshot = snapshot
        self.index = index

    def __repr__(self):
        return "Beta atom # {}".format(self.index)

    @property
    def xyz(self):
        return self.centroid

    @property
    def centroid(self):
        """
        Gets the centroid of this beta atom


        Returns
        -------

        centroid coordinate : np.ndarray
            shape = (3,)
        """
        return np.mean([alpha.centroid for alpha in self.alphas], axis=0)

    @property
    def alphas(self):
        for i in self.snapshot.beta_alpha_index_list[self.index]:
            yield Alpha(self.snapshot, i)

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
    def scores(self):
        return self.snapshot.beta_scores[self.index]

    @property
    def best_probe_type(self):
        return ['C', 'Br', 'F', 'Cl', 'I', 'OA', 'SA', 'N', 'P'][int(np.argmin(self.scores))]

    @property
    def score(self):
        """
        Get the score of this beta atom, which is the lowest vina score in all 9 probes.

        Returns
        -------
        float

        """
        return np.min(self.scores)

    @property
    def is_contact(self):
        return bool(self.snapshot.beta_contact[self.index])


class Pocket:

    def __init__(self, snapshot, index=None, beta_index=None):

        self.snapshot = snapshot
        self.void = False
        if index is not None:
            self.beta_index = self.snapshot.pocket_beta_index_list[index]
        elif beta_index is not None:
            self.beta_index = beta_index
        else:
            self.void = True


    @property
    def is_contact(self):
        if self.void:
            return False
        return np.any(self.snapshot.beta_contact[self.beta_index])

    @property
    def alphas(self):
        if self.void:
            return
        for i in self.beta_index:
            for j in self.snapshot.beta_alpha_index_list[i]:
                yield Alpha(snapshot=self.snapshot, index=j)

    @property
    def betas(self):
        if self.void:
            return
        for i in self.beta_index:
            yield Beta(snapshot=self.snapshot, index=i)

    @property
    def space(self):
        if self.void:
            return 0
        return np.sum([beta.space for beta in self.betas])

    @property
    def score(self):
        if self.void:
            return 0
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
        if self.void:
            return
        return np.mean([beta.centroid for beta in self.betas], axis=0)


class DPocket:

    def __init__(self, universe, beta_indices=None):
        """

        Parameters
        ----------
        universe: AS_Universe
        beta_indices
        """
        self.universe = universe

        self._betas = self.universe.map_beta(beta_indices) if beta_indices is not None else None

    def __len__(self):
        return self.n_pockets

    def __mul__(self, other):
        assert len(self.ss_vector) == len(other.ss_vector)

        return self.ss_vector * other.ss_vector

    def __add__(self, other):
        """

        Parameters
        ----------
        other: DPocket

        Returns
        -------
        combined_pocket: DPocket

        """
        assert isinstance(other, DPocket)

        combined_dpocket = DPocket(self.universe, beta_indices=None)
        combined_dpocket._betas = np.concatenate((self._betas, other._betas))

        return combined_dpocket

    @property
    def ss_vector(self):
        return np.array([len(self._betas[:, 1][self._betas[:, 0] == i]) / float(len(self._betas)) for i in
                         range(len(self.universe))])

    @property
    def variance(self):
        return np.sum((self.ss_vector - np.full(len(self.ss_vector), fill_value=1. / len(self.universe))) ** 2)

    @property
    def betas(self):
        for ss_idx, beta_idx in self._betas:
            yield Beta(self.universe[ss_idx], beta_idx)

    @property
    def pockets(self):
        """

        Returns
        -------
        Yields:
        Pocket


        """
        for i in range(len(self.universe)):
            beta_indices = self._betas[:, 1][self._betas[:, 0] == i]
            if len(beta_indices) > 0:
                yield Pocket(snapshot=self.universe[i], beta_index=beta_indices)
            else:
                yield Pocket(snapshot=self.universe[i])

    @property
    def scores(self):
        return np.array([pocket.score for pocket in self.pockets])

    @property
    def spaces(self):
        return np.array([pocket.space for pocket in self.pockets])

    @property
    def is_contact(self):
        return np.any([beta.is_contact for beta in self.betas])

    @property
    def centroid(self):
        return np.mean([beta.centroid for beta in self.betas], axis=0)

    @property
    def n_pockets(self):
        """
        Calculate how many pockets are in the d_pocket, which is equivalent to the number of snapshot the dpocket covers
        Returns
        -------
        int
        """
        return np.unique(self._betas[:, 0], return_counts=True).astype(int)
