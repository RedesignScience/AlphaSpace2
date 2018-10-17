import numpy as np
from .functions import _group, _getTetrahedronVolumes, _markInRange
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial import Delaunay, Voronoi
from .Cluster import _Pocket


class Snapshot:
    min_r = 3.2
    max_r = 5.4

    beta_cluster_dist = 1.3
    pocket_cluster_dist = 4.7

    contact_cutoff = 1.6

    def __init__(self):
        self.residue_names = None
        self.elements = None
        self.atom_names = None

        self._alpha_xyz = None
        self._alpha_space = None
        self._alpha_contact = None

        self._beta_xyz = None
        self._beta_space = None
        self._beta_contact = None
        self._beta_scores = None

        self._beta_dpocket_index = None

        self._pocket_xyz = None
        self._pocket_space = None
        self._pocket_contact = None

        self._alpha_lining = None
        self._alpha_radii = None

        self._beta_alpha_index_list = None
        self._pocket_beta_index_list = None

    def genAlphas(self, receptor):
        """


        Parameters
        ----------
        receptor : mdtraj.trajectory
        """
        self.residue_names, self.elements, self.atom_names = zip(
            *[(atom.residue.name, atom.element.symbol, atom.name) for atom in receptor.top.atoms])

        protein_coords = receptor.xyz[0]

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

        self._alpha_radii = np.take(raw_alpha_sphere_radii, filtered_alpha_idx)

        self._alpha_lining = np.take(raw_alpha_lining_idx, filtered_alpha_idx, axis=0)

        alpha_lining_xyz = np.take(protein_coords, self._alpha_lining, axis=0).astype(np.float32)

        self._alpha_space = _getTetrahedronVolumes(alpha_lining_xyz)

        self._alpha_xyz = np.take(raw_alpha_xyz, filtered_alpha_idx, axis=0)

    def genBetas(self):
        assert self._alpha_xyz is not None
        zmat = linkage(self._alpha_xyz, method='average')

        alpha_beta_label = fcluster(zmat, self.beta_cluster_dist / 10, criterion='distance') - 1

        self._beta_alpha_index_list = _group(alpha_beta_label)

        self._beta_xyz = [None] * (max(alpha_beta_label) + 1)
        self._beta_space = [None] * (max(alpha_beta_label) + 1)
        for i, indices in enumerate(self._beta_alpha_index_list):
            self._beta_xyz[i] = np.mean(self._alpha_xyz[indices], axis=0)
            self._beta_space[i] = np.sum(self._alpha_space[indices], axis=0)

        self._beta_xyz = np.array(self._beta_xyz)
        self._beta_space = np.array(self._beta_space)

    def genPockets(self):
        zmat = linkage(self._beta_xyz, method='average')
        beta_pocket_label = fcluster(zmat, self.pocket_cluster_dist / 10, criterion='distance') - 1

        self._pocket_beta_index_list = _group(beta_pocket_label)
        self._pocket_xyz = [None] * (max(beta_pocket_label) + 1)
        self._pocket_space = [None] * (max(beta_pocket_label) + 1)
        for i, indices in enumerate(self._pocket_beta_index_list):
            self._pocket_xyz[i] = np.mean(self._beta_xyz[indices], axis=0)
            self._pocket_space[i] = np.sum(self._beta_space[indices], axis=0)

    def genBScore(self, receptor):
        from .VinaScoring import _pre_process_pdbqt, _get_probe_score

        if hasattr(receptor, 'adv_atom_types'):

            prot_types, hp_type, acc_type, don_type = _pre_process_pdbqt(receptor)

            self._beta_scores = _get_probe_score(probe_coords=self._beta_xyz * 10, prot_coord=receptor.xyz[0] * 10,
                                                 prot_types=prot_types,
                                                 hp_type=hp_type,
                                                 acc_type=acc_type, don_type=don_type)
            print("Vina Atom Type found, calculating BScore")
        else:
            self._beta_scores = np.zeros(len(self._beta_xyz), dtype=np.float)
            print("No Vina Atom Type found")

    def calculateContact(self, coords):
        """
        Mark alpha/beta/pocket atoms as contact with in cutoff of ref points.

        _Beta atom and pocket atoms are counted as contact if any of their child alpha atoms is in contact.

        Parameters
        ----------
        coords: np.array shape = (n,3)

        Returns
        -------
        """

        self._alpha_contact = _markInRange(self._alpha_xyz, ref_points=coords, cutoff=self.contact_cutoff / 10)
        self._beta_contact = self._beta_contact = np.array(
            [np.any(self._alpha_contact[alpha_indices]) for alpha_indices in self._beta_alpha_index_list])
        self._pocket_contact = np.array(
            [np.any(self._beta_contact[beta_indices]) for beta_indices in self._pocket_beta_index_list])

    def run(self, receptor, binder=None):

        self.genAlphas(receptor)
        self.genBetas()
        self.genPockets()

        self.genBScore(receptor)

        if binder is not None:
            self.calculateContact(coords=binder.xyz[0])

    @property
    def pockets(self):
        for i in range(len(self._pocket_xyz)):
            yield _Pocket(self, i)

    @property
    def betas(self):
        for p in self.pockets:
            for b in p.betas:
                yield b

    @property
    def alphas(self):
        for p in self.pockets:
            for a in p.alphas:
                yield a
