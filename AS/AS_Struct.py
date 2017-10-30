import numpy as np
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform,jaccard,pdist
from AS_Cluster import AS_D_Pocket,AS_Cluster
from AS_Funct import _tessellation




class AS_Structure:
    def __init__(self, trajectory: object, structure_type: int = 2, parent: object = None) -> object:

        """
        Container for structure trajectory and topology in a AS_Session
        :param trajectory: MDtraj trajectory object
        :param structure_type: int 0 for receptor, 1 for binder, 2 for other
        """

        # 0 for receptor, 1 for binder, 2 for unassigned
        self.structure_type = structure_type
        self._clusters = {}
        self.trajectory = trajectory
        self.parent = parent
        self.universe = parent
        self.config = self.parent.config
        self.contact_cluster = [[None for i in range(self.n_residues)] for j in range(self.n_frames)]

    def __repr__(self):
        return "{} Structure with {} frames, {} residues, {} atoms".format(
            ['Receptor','Binder','Misc.'][self.structure_type],self.n_frames,self.n_residues,self.n_atoms)

    @property
    def is_polar(self):
        return np.array(
            [(str(atom.element) in ['nitrogen', 'oxygen', 'sulfur']) for atom in self.topology._atoms])

    @property
    def top(self):
        return self.trajectory.topology

    @property
    def traj(self):
        return self.trajectory

    @property
    def n_atoms(self):
        return self.trajectory.n_atoms

    @property
    def n_frames(self):
        return self.trajectory.n_frames

    @property
    def n_residues(self):
        return self.topology.n_residues

    @property
    def topology(self):
        return self.trajectory.topology


    def cluster(self, snapshot_idx):
        return self._clusters[snapshot_idx]

    @property
    def clusters(self) -> object:
        for key in sorted(self._clusters.keys()):
            yield self._clusters[key]

    @property
    def n_clusters(self) -> int:
        return self.trajectory.n_frames

    @property
    def residues(self):
        """
        Residue iterator
        :return: iter
        """
        return self.topology.residues

    @property
    def atoms(self):
        """
        Atom iterator
        :return: iter
        """
        for atom in self.top._atoms:
            yield atom


    @property
    def __len__(self):
        """
        Returns number of frames
        :return: int
        """
        return self.n_frames

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

    def calculate_contact(self, binder, snapshot_idx=0):
        """
        Calculate the contact index of the alpha cluster against the designated binder.
        The contact distance cutoff can be set in config
        :param binder: object, AS_Struct
        :param snapshot_idx: int
        """
        self._clusters[snapshot_idx]._get_contact_space()

    def generate_cluster(self, snapshot_idx=0):
        """
        Perform tessellation of a receptor snapshot
        :param snapshot_idx: int
        """
        self._clusters[snapshot_idx] = AS_Cluster(self, snapshot_idx)
        return self._clusters[snapshot_idx]

    def assign_binder_contact_pocket(self, AS_Cluster, snapshot_idx):
        contact_matrix = AS_Cluster._get_contact_list(self.trajectory[snapshot_idx])
        for residue in self.topology.residues:
            self.contact_cluster[snapshot_idx][residue.index] = AS_Cluster._get_pockets_by_binder_contact(
                contact_matrix=contact_matrix,
                binder_residue=residue)

    def get_residue_contact_pocket(self, AS_Cluster, residue_index):
        if self.contact_cluster[AS_Cluster.snapshot_idx][0] is None:
            self.assign_binder_contact_pocket(AS_Cluster, AS_Cluster.snapshot_idx)
        contact_pocket_index = self.contact_cluster[AS_Cluster.snapshot_idx][residue_index]
        for i in contact_pocket_index:
            yield AS_Cluster.pocket(i)

    def _gen_d_pockets(self) -> object:
        """
        Generate d-pocket dictionary of list of indices
        :return: dict of d_pockets
        """
        if len(self._clusters) != self.n_frames:
            raise Exception("All Frames must be processed before d pocket generation")

        pockets = []
        for i in range(self.n_frames):
            pockets.extend(list(self.cluster(snapshot_idx=i).pockets))

        lining_atoms = [p.get_lining_atoms() for p in pockets]
        lining_dist = np.zeros((len(lining_atoms),len(lining_atoms)))

        for i in range(len(lining_atoms)):
            for j in range(i,len(lining_atoms)):
                lining_dist[i,j] = lining_dist[j,i] = len(lining_atoms[i].symmetric_difference(lining_atoms[j]))/len(lining_atoms[i].union(lining_atoms[j]))
        lining_pdist = squareform(lining_dist)

        clustered_list = list(fcluster(Z=linkage(lining_pdist, method='complete'), t=self.config.dpocket_cluster_cutoff,
                                       criterion='distance'))

        d_pockets_idx_dict = {}
        for pocket, d_pocket_idx in zip(pockets,clustered_list):
            snapshot_idx = pocket.cluster.snapshot_idx
            pocket_idx = pocket.index
            if d_pocket_idx not in d_pockets_idx_dict:
                d_pockets_idx_dict[d_pocket_idx] = []
            d_pockets_idx_dict[d_pocket_idx].append((snapshot_idx,pocket_idx))

        for idx,item in d_pockets_idx_dict.items():
            yield AS_D_Pocket(parent_universe=self.parent,pocket_indices=item)
