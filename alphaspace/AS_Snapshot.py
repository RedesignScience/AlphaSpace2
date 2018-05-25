import numpy as np
from .AS_Funct import group, getTetrahedronVolumes, mark_in_range
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial import Delaunay, Voronoi
from .AS_Cluster import Alpha, Beta, Pocket


class AS_Snapshot:
    min_r = 3.2
    max_r = 5.4

    beta_cluster_dist = 1.6

    pocket_cluster_dist = 4.7

    def __init__(self, traj=None, snapshot_idx=0):
        self.residue_names = None
        self.elements = None
        self.atom_names = None

        self.alpha_xyz = None
        self.alpha_space = None
        self._alpha_contact = None

        self.beta_xyz = None
        self.beta_space = None
        self._beta_contact = None

        self.pocket_xyz = None
        self.pocket_space = None
        self._pocket_contact = None

        self.alpha_lining = None
        self.pocket_xyz = None
        self.alpha_radii = None

        self.beta_alpha_index_list = None
        self.pocket_beta_index_list = None

        self.beta_scores = None

        self.snapshot_idx = snapshot_idx

        if traj is not None:
            self.tessellation(traj, snapshot_idx)

    @property
    def n_betas(self):
        return len(self.beta_xyz)


    def tessellation(self, traj, snapshot_idx=0):
        """


        Parameters
        ----------
        traj : mdtraj.trajectory

        Returns
        -------
        alpha_xyz
        alpha _lining
        alpha_radii
        """
        self.residue_names = [atom.residue.name for atom in traj.top.atoms]

        self.elements = [atom.element.symbol for atom in traj.top.atoms]

        self.atom_names = [atom.name for atom in traj.top.atoms]

        protein_coords = traj.xyz[snapshot_idx]

        raw_alpha_lining_idx = Delaunay(protein_coords).simplices

        # Take coordinates from xyz file
        raw_alpha_lining_xyz = np.take(protein_coords, raw_alpha_lining_idx[:, 0].flatten(), axis=0)

        # generate alpha atom coordinates
        raw_alpha_xyz = Voronoi(protein_coords).vertices
        # Calculate alpha sphere radii
        raw_alpha_sphere_radii = np.linalg.norm(raw_alpha_lining_xyz - raw_alpha_xyz, axis=1)

        # Filter the data based on radii cutoff
        filtered_alpha_idx = np.where(np.logical_and(self.min_r / 10.0 <= raw_alpha_sphere_radii,
                                                     raw_alpha_sphere_radii <= self.max_r / 10.0))[0]

        self.alpha_radii = np.take(raw_alpha_sphere_radii, filtered_alpha_idx)

        self.alpha_lining = np.take(raw_alpha_lining_idx, filtered_alpha_idx, axis=0)

        alpha_lining_xyz = np.take(protein_coords, self.alpha_lining, axis=0).astype(np.float32)

        self.alpha_space = getTetrahedronVolumes(alpha_lining_xyz)

        self.alpha_xyz = np.take(raw_alpha_xyz, filtered_alpha_idx, axis=0)

        return self.alpha_xyz, self.alpha_lining, self.alpha_radii

    def _gen_beta(self, dist=None):
        zmat = linkage(self.alpha_xyz, method='average')

        dist = self.beta_cluster_dist if dist is None else dist

        alpha_beta_label = fcluster(zmat, dist / 10, criterion='distance') - 1

        self.beta_alpha_index_list = group(alpha_beta_label)

        self.beta_xyz = [None] * (max(alpha_beta_label) + 1)
        self.beta_space = [None] * (max(alpha_beta_label) + 1)
        for i, indices in enumerate(self.beta_alpha_index_list):
            self.beta_xyz[i] = np.mean(self.alpha_xyz[indices], axis=0)
            self.beta_space[i] = np.sum(self.alpha_space[indices], axis=0)

        self.beta_xyz = np.array(self.beta_xyz)
        self.beta_space = np.array(self.beta_space)

    def _gen_pocket(self, dist=None):
        zmat = linkage(self.beta_xyz, method='average')

        dist = self.pocket_cluster_dist if dist is None else dist

        beta_pocket_label = fcluster(zmat, dist / 10, criterion='distance') - 1

        self.pocket_beta_index_list = group(beta_pocket_label)

        self.pocket_xyz = [None] * (max(beta_pocket_label) + 1)
        self.pocket_space = [None] * (max(beta_pocket_label) + 1)
        for i, indices in enumerate(self.pocket_beta_index_list):
            self.pocket_xyz[i] = np.mean(self.beta_xyz[indices], axis=0)
            self.pocket_space[i] = np.sum(self.beta_space[indices], axis=0)

    def screen_by_contact(self, ref_points, cutoff=3.6):
        """
        Mark alpha/beta/pocket atoms as contact with in cutoff of ref points.

        Beta atom and pocket atoms are counted as contact if any of their child alpha atoms is in contact.

        Parameters
        ----------
        ref_points: np.array shape = (n,3)
        cutoff: float

        Returns
        -------
        """

        self._alpha_contact = mark_in_range(self.alpha_xyz, ref_points=ref_points, cutoff=cutoff / 10)
        self._pocket_contact = None
        self._beta_contact = None

    @property
    def alpha_contact(self):
        if self._alpha_contact is None:
            raise Exception("No contact has been set yet, use screen_by_contact method")
        return self.alpha_contact

    @property
    def beta_contact(self):
        if self._beta_contact is None:
            self._beta_contact = np.array(
                [np.any(self._alpha_contact[alpha_indices]) for alpha_indices in self.beta_alpha_index_list])
        return self._beta_contact

    @property
    def pocket_contact(self):
        if self._pocket_contact is None:
            self._pocket_contact = np.array(
                [np.any(self.beta_contact[beta_indices]) for beta_indices in self.pocket_beta_index_list])
        return self._pocket_contact

    @property
    def pockets(self):
        for i in range(len(self.pocket_xyz)):
            yield Pocket(self, i)
