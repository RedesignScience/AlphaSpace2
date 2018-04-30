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
        return self.snapshot.alpha_space[self.index] * 100

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

    def __init__(self, snapshot, index):

        self.snapshot = snapshot
        self.index = index

    @property
    def is_contact(self):
        return self.snapshot.pocket_contact[self.index]

    @property
    def alphas(self):
        for i in self.snapshot.pocket_beta_index_list[self.index]:
            for j in self.snapshot.beta_alpha_index_list[i]:
                yield Alpha(snapshot=self.snapshot, index=j)
    @property
    def betas(self):

        for i in self.snapshot.pocket_beta_index_list[self.index]:
            yield Beta(snapshot=self.snapshot, index=i)


    @property
    def space(self):
        return np.sum([beta.space for beta in self.betas])

    @property
    def score(self):
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

        return np.mean([beta.centroid for beta in self.betas], axis=0)
