"""
This file contains the container classes for cluster based objects.

AS_Data: A numpy.ndarray inheritance that stores and enumerates data

AS_Pocket: A mask container for pocket information, gets data from AS_Data

AS_AlphaAtom: A mast container for alpha atom, gets data from AS_Data

AS_BetaATOM: same as alpha atom

"""


from itertools import chain

import numpy as np


class _Alpha:
    """
    This is the alpha atom container, which is the constituents of the pocket
    object. This object only contains reference of alpha atom id, and all info
    is stored in the array you can access by methods.
    """

    def __init__(self, snapshot, index):
        """

        Parameters
        ----------
        snapshot: alpsphace.Snapshot
        index: int
        """

        self.snapshot = snapshot
        self.index = index

    def __repr__(self):
        return "AlphaAtom {}".format(self.index)

    def __sub__(self, other):
        """
        Use subtraction to get the distance between two alpha atoms

        Parameters
        ----------
        other : _Alpha

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
        return self.snapshot._alpha_xyz[self.index]

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
        return self.snapshot._alpha_space[self.index]

    @property
    def nonpolar_space(self):
        return self.snapshot._alpha_space[self.index] * self.snapshot._alpha_space_nonpolar_ratio[self.index]

    @property
    def lining_atoms_idx(self):
        """
        Get the atom index of all the lining atoms. Each alpha atom has four lining atoms.


        Returns
        -------

        indices : np.ndarray
            size (4)

        """
        return self.snapshot._alpha_lining[self.index]

    @property
    def isContact(self) -> bool:
        """
        Check if it's in contact with any binder atoms

        Returns
        -------

        bool

        """
        return bool(self.snapshot._alpha_contact[self.index])


class _Beta:
    """
    This is the container for _Beta atom,
    which is simply a collection of alpha atoms.

    It belongs to the AS_Pocket object.
    """

    def __init__(self, snapshot, index):
        self.snapshot = snapshot
        self.index = index

    def __repr__(self):
        return "_Beta atom # {}".format(self.index)

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
        for i in self.snapshot._beta_alpha_index_list[self.index]:
            yield _Alpha(self.snapshot, i)

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
    def nonpolar_space(self):
        return np.sum([alpha.nonpolar_space for alpha in self.alphas])

    @property
    def scores(self):
        return self.snapshot._beta_scores[self.index]

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
    def isContact(self):
        return bool(self.snapshot._beta_contact[self.index])

    @property
    def occupiedSpace(self):
        return np.sum(
            np.array([alpha.space for alpha in self.alphas]) *
            np.array([alpha.isContact for alpha in self.alphas])
        )

    @property
    def occupiedNonpolarSpace(self):
        return np.sum(
            np.array([alpha.nonpolar_space for alpha in self.alphas]) *
            np.array([alpha.isContact for alpha in self.alphas])
        )

    @property
    def occupancy(self):
        return self.occupiedSpace / self.space

    @property
    def occupancy_nonpolar(self):
        return self.occupiedNonpolarSpace / self.nonpolar_space


class _Pocket:
    def __init__(self, snapshot, pocketIndex=None, alphaIndex=None, betaIndex=None):
        self.snapshot = snapshot
        self.isEmpty = False
        self.index = pocketIndex
        if betaIndex is not None:
            self.beta_index = betaIndex
            self.alpha_index = list(
                chain(*[self.snapshot._beta_alpha_index_list[i] for i in betaIndex]))
        elif alphaIndex is not None:
            self.alpha_index = alphaIndex
        elif pocketIndex is not None:
            self.alpha_index = self.snapshot._pocket_alpha_index_list[pocketIndex]
        else:
            self.alpha_index = []
            self.isEmpty = True

    def __repr__(self):
        return "Pocket {}".format(self.index)

    @property
    def isContact(self):
        if self.isEmpty:
            return False
        if self.index is None:
            return np.any(self.snapshot._alpha_contact[self.alpha_index])
        else:
            return bool(self.snapshot._pocket_contact[self.index])

    @property
    def occupiedSpace(self):
        if self.isEmpty:
            return 0.0
        return np.sum(
            np.array([alpha.space for alpha in self.alphas]) * np.array([alpha.isContact for alpha in self.alphas]))

    @property
    def occupiedNonpolarSpace(self):
        if self.isEmpty:
            return 0.0
        return np.sum(
            np.array([alpha.nonpolar_space for alpha in self.alphas]) * np.array([alpha.isContact for alpha in self.alphas]))

    @property
    def occupancy(self):
        return self.occupiedSpace / self.space

    @property
    def occupancy_nonpolar(self):
        return self.occupiedNonpolarSpace / self.nonpolar_space

    @property
    def alphas(self):
        if self.isEmpty:
            return
        for i in self.alpha_index:
            yield _Alpha(snapshot=self.snapshot, index=i)

    @property
    def betas(self):
        if self.isEmpty:
            return
        elif hasattr(self, 'beta_index'):
            for i in self.beta_index:
                yield _Beta(snapshot=self.snapshot, index=i)
        else:
            for i in self.snapshot._pocket_beta_index_list[self.index]:
                yield _Beta(snapshot=self.snapshot, index=i)

    @property
    def space(self):
        if self.isEmpty:
            return 0
        return np.sum([alpha.space for alpha in self.alphas])

    @property
    def nonpolar_space(self):
        if self.isEmpty:
            return 0
        return np.sum([alpha.nonpolar_space for alpha in self.alphas])


    @property
    def score(self):
        if self.isEmpty:
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
        if self.isEmpty:
            return
        return np.mean([alphas.centroid for alphas in self.alphas], axis=0)

    @property
    def lining_atoms_idx(self):
        return np.unique([alpha.lining_atoms_idx for alpha in self.alphas])


class _DPocket:

    def __init__(self, trajectory, beta_indices=None):
        self.trajectory = trajectory
        self._betas = self.trajectory.map_beta(
            beta_indices) if beta_indices is not None else None

    def __len__(self):
        return self.n_pockets

    def __add__(self, other):
        """

        Parameters
        ----------
        other: _DPocket

        Returns
        -------
        combined_pocket: _DPocket

        """
        assert isinstance(other, _DPocket)

        combined_dpocket = _DPocket(self.trajectory, beta_indices=None)
        combined_dpocket._betas = np.concatenate((self._betas, other._betas))

        return combined_dpocket

    @property
    def isContact(self):
        for pocket in self.pockets:
            if pocket.isContact:
                return True
        return False

    def __repr__(self):
        return "DPocket from {} snapshots {} beta atoms".format(str(len(self.trajectory)), str(len(self._betas)))

    @property
    def betas(self):
        for ss_idx, beta_idx in self._betas:
            yield _Beta(self.trajectory[ss_idx], beta_idx)

    @property
    def pockets(self):
        """

        Returns
        -------
        Yields:
        _Pocket
        """
        for i in range(len(self.trajectory)):
            beta_indices = np.array(self._betas[:, 1][self._betas[:, 0] == i])
            if len(beta_indices) > 0:
                yield _Pocket(snapshot=self.trajectory[i], betaIndex=beta_indices)
            else:
                yield _Pocket(snapshot=self.trajectory[i])

    @property
    def scores(self):
        return np.array([pocket.score for pocket in self.pockets])

    @property
    def spaces(self):
        return np.array([pocket.space for pocket in self.pockets])

    @property
    def nonpolar_spaces(self):
        return np.array([pocket.nonpolar_space for pocket in self.pockets])

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

