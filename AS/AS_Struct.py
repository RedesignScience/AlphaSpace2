from AS.AS_Cluster import *


class AS_Structure:
    def __init__(self, trajectory, structure_type=2, parent=None):

        """
        Container for structure trajectory and topology in a AS_Session
        :param trajectory: MDtraj trajectory object
        :param structure_type: int 0 for receptor, 1 for binder, 2 for other
        """

        # 0 for receptor, 1 for binder, 3 for unassigned
        self.structure_type = structure_type
        self._clusters = {}
        self.trajectory = trajectory
        self.parent = parent
        self.config = self.parent.config
        self.is_polar = np.array(
            [(str(atom.element) in ['nitrogen', 'oxygen', 'sulfur']) for atom in self.topology.atoms])
        self.contact_cluster = [[None for i in range(self.n_residues)] for j in range(self.n_frames)]

    def __len__(self):
        """
        Returns number of frames
        :return: int
        """
        return self.n_frames

    def __repr__(self):
        if self.structure_type == 0:
            return "Receptor part trajectory of {} residues in {} snapshots".format(self.topology.n_residues, len(self))
        elif self.structure_type == 1:
            return "Binder part trajectory of {} residues in {} snapshots".format(self.topology.n_residues, len(self))
        else:
            return "Unknown part trajectory of {} residues in {} snapshots".format(self.topology.n_residues, len(self))

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
        return self.clusters[snapshot_idx]

    @property
    def clusters(self):
        return self._clusters

    @property
    def n_clusters(self):
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
        return self.topology.atoms

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

    def __bool__(self):
        """
        Check if empty
        :return: bool
        """
        return len(self) > 0

    def calculate_contact(self, binder, snapshot_idx=0):
        """
        Calculate the contact index of the alpha cluster against the designated binder.
        The contact distance cutoff can be set in config
        :param binder: object, AS_Struct
        :param snapshot_idx: int
        """
        if binder:
            self.clusters[snapshot_idx]._get_contact_space()

    def generate_cluster(self, snapshot_idx=0):
        """
        Perform tessellation of a receptor snapshot
        :param snapshot_idx: int
        """

        self.clusters[snapshot_idx] = AS_Cluster(self, snapshot_idx)

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
