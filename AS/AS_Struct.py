from AS.AS_Cluster import *


class AS_Structure:
    def __init__(self,trajectory,structure_type=2):

        """
        Container for structure trajectory and topology in a AS_Session
        :param trajectory: MDtraj trajectory object
        :param structure_type: int 0 for receptor, 1 for binder, 2 for other
        """

        # 0 for receptor, 1 for binder, 3 for unassigned
        self.structure_type = structure_type
        self.topology = trajectory.topology
        self.trajectory = trajectory
        self.n_atoms = self.topology.n_atoms
        self.n_frames = self.trajectory.n_frames
        self.n_residues = self.topology.n_residues
        self.clusters = [None for _ in range(len(self))]

        self.is_polar = np.array(
            [(str(atom.element) in ['nitrogen', 'oxygen', 'sulfur']) for atom in self.topology.atoms])

        self.contact_cluster = [[None for i in range(self.n_residues) ] for j in range(self.n_frames)]


    def __len__(self):
        """
        Returns number of frames
        :return: int
        """
        return self.n_frames

    def get_cluster(self,snapshot_idx):
        return self.clusters[snapshot_idx]

    def __getitem__(self, item):
        return self.get_cluster(item)

    def __repr__(self):
        if self.structure_type == 0:
            return "Receptor part trajectory of {} residues in {} snapshots".format(self.topology.n_residues,len(self))
        elif self.structure_type == 1:
            return "Binder part trajectory of {} residues in {} snapshots".format(self.topology.n_residues,len(self))
        else:
            return "Unknown part trajectory of {} residues in {} snapshots".format(self.topology.n_residues,len(self))

    def residues(self):
        """
        Residue iterator
        :return: iter
        """
        return self.topology.residues()

    def atoms(self):
        """
        Atom iterator
        :return: iter
        """
        return self.topology.atoms()

    def residue(self,idx):
        """
        Gives a residue with idx
        :param idx: int
        :return: Residue
        """
        return self.topology.residue(idx)

    def atom(self,idx):
        """
        Gives an atom with idx
        :param idx: int
        :return: Atom
        """
        return self.topology.atom(idx)

    def __bool__(self):
        """
        Check if empty
        :return: bool
        """
        if len(self) > 0:
            return True
        else:
            return False

    def calculate_contact(self,binder,snapshot_idx=0):
        """
        Calculate the contact index of the alpha cluster against the designated binder.
        The contact distance cutoff can be set in config
        :param binder: object, AS_Struct
        :param snapshot_idx: int
        """
        if binder:
            self.clusters[snapshot_idx].calculate_contact(binder.trajectory.xyz[snapshot_idx])


    def assign_cluster_contact(self,AS_Cluster,snapshot_idx):
        contact_matrix = AS_Cluster._get_binder_contact(self.trajectory[snapshot_idx])
        for residue in self.topology.residues:
            self.contact_cluster[snapshot_idx][residue.index] = AS_Cluster._take_contact_list(binder_pocket_contact_matrix=contact_matrix,
                                                                                              binder_residue=residue)


    def get_residue_contact_pocket(self,AS_Cluster,residue_index):
        if self.contact_cluster[AS_Cluster.snapshot_idx][0] is None:
            self.assign_cluster_contact(AS_Cluster,AS_Cluster.snapshot_idx)
        contact_pocket_index = self.contact_cluster[AS_Cluster.snapshot_idx][residue_index]
        for i in contact_pocket_index:
            yield AS_Cluster.pocket(i)








