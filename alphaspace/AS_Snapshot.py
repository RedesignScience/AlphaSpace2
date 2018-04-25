from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial import Delaunay, Voronoi
import numpy as np
from alphaspace.AS_Funct import group, getTetrahedronVolumes, find_in_range


class AS_Snapshot:
    min_r = 3.2
    max_r = 5.4

    beta_cluster_dist = 1.8

    pocket_cluster_dist = 4.7

    def __init__(self, traj=None, snapshot_idx=0):
        self.residue_names = None
        self.elements = None
        self.atom_names = None

        self.alpha_xyz = None
        self.alpha_space = None
        self.contact_alpha = []

        self.beta_xyz = None
        self.beta_space = None
        self.contact_beta = []

        self.pocket_xyz = None
        self.pocket_space = None
        self.contact_pocket = []

        self.alpha_lining = None
        self.alpha_xyz = None
        self.alpha_radii = None

        self.beta_alpha_index_list = None

        self.beta_scores = None

        if traj is not None:
            self.tessellation(traj, snapshot_idx)

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
        self.residue_names = [atom.element.number for atom in traj.top.atoms]

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
