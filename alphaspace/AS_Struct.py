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
        self.config = self.parent.config
        # self.contact_cluster = [[None for i in range(self.n_residues)] for j in range(self.n_frames)]
        self._data = None
        self._pockets = None

    def _tessellation(self, config, snapshot_idx: int) -> np.ndarray:
        # Generate Raw Tessellation simplexes
        raw_alpha_lining_idx = Delaunay(self.trajectory.xyz[snapshot_idx]).simplices
        # Take coordinates from xyz file
        raw_alpha_lining_xyz = np.take(self.trajectory.xyz[snapshot_idx],raw_alpha_lining_idx[:,0].flatten(),axis=0)

        # generate alpha atom coordinates
        raw_alpha_xyz = Voronoi(self.traj.xyz[snapshot_idx]).vertices

        # Calculate alpha sphere radii
        raw_alpha_sphere_radii = np.linalg.norm(raw_alpha_lining_xyz - raw_alpha_xyz,axis=1)

        # Filter the data based on radii cutoff
        filtered_alpha_idx = np.where(np.logical_and(config.min_r / 10.0 <= raw_alpha_sphere_radii,
                                                     raw_alpha_sphere_radii <= config.max_r / 10.0))[0]

        filtered_alpha_radii = np.take(raw_alpha_sphere_radii,filtered_alpha_idx)

        alpha_lining = np.take(raw_alpha_lining_idx,filtered_alpha_idx,axis=0)

        filtered_alpha_xyz = np.take(raw_alpha_xyz,filtered_alpha_idx,axis=0)

        # cluster the remaining vertices to assign index of belonging pockets
        zmat = linkage(filtered_alpha_xyz,method='average')

        alpha_pocket_index = fcluster(zmat,self.config.clust_dist / 10,
                                      criterion='distance') - 1  # because cluster index start from 1

        # Load trajectories
        filtered_lining_xyz = np.take(self.traj.xyz[snapshot_idx],alpha_lining,axis=0)
        # calculate the polarity of alpha atoms
        _total_score = np.array(
                [getTetrahedronVolume(i) for i in
                 filtered_lining_xyz]) * 1000  # here the 1000 is to convert nm^3 to A^3

        atom_sasa = getSASA(self.trajectory[snapshot_idx])
        atom_covered_sasa = getSASA(self.traj[snapshot_idx], filtered_alpha_xyz)
        pocket_sasa = np.take(atom_sasa - atom_covered_sasa,alpha_lining)

        is_polar = np.array(
                [(str(atom.element) in ['nitrogen','oxygen','sulfur']) for atom in self.topology._atoms])
        polar_ratio = np.average(np.take(is_polar,alpha_lining),axis=1,weights=pocket_sasa)

        _polar_score = _total_score * polar_ratio

        _nonpolar_score = _total_score - _polar_score

        # for item in (np.zeros((len(alpha_pocket_index),1)),
        #                        np.full((len(alpha_pocket_index),1),snapshot_idx),
        #                        filtered_alpha_xyz,
        #                        alpha_lining,
        #                        np.expand_dims(_polar_score,axis=1),
        #                        np.expand_dims(_nonpolar_score, axis=1),
        #                        np.ones((len(alpha_pocket_index),1)),
        #                        np.zeros((len(alpha_pocket_index),1)),
        #                        np.expand_dims(alpha_pocket_index, axis=1),
        #                        np.expand_dims(filtered_alpha_radii,axis=1)
        #                        ):
        #     print(item.shape)

        data = np.concatenate((np.zeros((len(alpha_pocket_index),1)),
                               np.full((len(alpha_pocket_index),1),snapshot_idx),
                               filtered_alpha_xyz,
                               alpha_lining,
                               np.expand_dims(_polar_score,axis=1),
                               np.expand_dims(_nonpolar_score, axis=1),
                               np.ones((len(alpha_pocket_index),1)),
                               np.zeros((len(alpha_pocket_index),1)),
                               np.expand_dims(alpha_pocket_index, axis=1),
                               np.expand_dims(filtered_alpha_radii,axis=1)
                               ),axis=-1)

        """
        0       idx
        1       snapshot_idx
        2       x
        3       y
        4       z
        5       lining_atom_idx_1
        6       lining_atom_idx_1
        7       lining_atom_idx_1
        8       lining_atom_idx_1
        9       polar_score 0
        10      nonpolar_score 0
        11      is_active 1
        12      is_contact 0
        13      pocket_idx
        14      radii
        """
        return data

    @property
    def is_polar(self):
        return np.array([(str(atom.element) in ['nitrogen','oxygen','sulfur']) for atom in self.topology._atoms])


    def _combine_data(self,data_list):
        assert type(data_list) == list
        assert type(data_list[0]) == np.ndarray
        data_list.sort(key = lambda d:d[0,1])
        data = np.concatenate(data_list)
        data[:, 0] = np.arange(0, len(data), dtype=int)
        self._data = AS_Data(data, self)


    def _gen_pockets(self):
        assert self._data is not None
        self._pockets = {i: {} for i in range(self.n_frames)}
        pocket_snapshot_dict = self._data[:,[0,1,13]]
        for i in range(self.n_frames):
            reversed_dict = defaultdict(list)
            current_frame_pocket_idx = pocket_snapshot_dict[pocket_snapshot_dict[:,1] == i]
            for idx,ss_idx,p_idx in current_frame_pocket_idx:
                reversed_dict[p_idx].append(idx)
            self._pockets[i] = reversed_dict

    def pockets(self,snapshot_idx=0):
        if self._pockets is None:
            self._gen_pockets()
        for pocket_idx,pocket_content in self._pockets[snapshot_idx].items():
            yield AS_Pocket(pocket_content,snapshot_idx,pocket_idx,self)



    def __repr__(self):
        return "{} Structure with {} frames, {} residues, {} atoms".format(
                ['Receptor','Binder','Misc.'][self.structure_type],self.n_frames,self.n_residues,self.n_atoms)

    @property
    def is_polar(self):
        return np.array(
                [(str(atom.element) in ['nitrogen','oxygen','sulfur']) for atom in self.topology._atoms])

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
    def n_snapshots(self):
        return self.trajectory.n_frames

    @property
    def n_residues(self):
        return self.topology.n_residues

    @property
    def topology(self):
        return self.trajectory.topology

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
        :return: object atom
        """
        return self.topology.atom(idx)


    def n_alphas(self,snapshot_idx = 0,active_only = False):
        """
        return the number of alpha atoms in a particular snapshot
        :param snapshot_idx: int
        :param active_only: bool
        :return: int
        """
        alpha_idx = self._data.snapshot_alpha_idx(snapshot_idx=snapshot_idx)
        if active_only:
            return len(np.where(self._data.is_active(alpha_idx)))
        else:
            return len(alpha_idx)




    def calculate_contact(self,snapshot_idx=0):
        """
        Calculate the contact index of the alpha cluster against the designated binder.
        The contact distance cutoff can be set in config
        :param snapshot_idx: int
        """
        alpha_idx = self._data.snapshot_alpha_idx(snapshot_idx)
        snapshot_cluster_coord_matrix = self._data.xyz(alpha_idx)
        binder_coords = self.universe.binder.trajectory.xyz[snapshot_idx]
        contact_alpha = getIfContact(snapshot_cluster_coord_matrix, binder_coords, self.config.hit_dist)[0]
        contact_alpha_idx = alpha_idx[contact_alpha]
        self._data[contact_alpha_idx,12] = 1

    # def assign_binder_contact_pocket(self,AS_Cluster,snapshot_idx):
    #     contact_matrix = AS_Cluster._get_contact_list(self.trajectory[snapshot_idx])
    #     for residue in self.topology.residues:
    #         self.contact_cluster[snapshot_idx][residue.index] = AS_Cluster._get_pockets_by_binder_contact(
    #                 contact_matrix=contact_matrix,
    #                 binder_residue=residue)

    # def get_residue_contact_pocket(self,AS_Cluster,residue_index):
    #     if self.contact_cluster[AS_Cluster.alpha_snapshot_idx][0] is None:
    #         self.assign_binder_contact_pocket(AS_Cluster, AS_Cluster.alpha_snapshot_idx)
    #     contact_pocket_index = self.contact_cluster[AS_Cluster.alpha_snapshot_idx][residue_index]
    #     for i in contact_pocket_index:
    #         yield AS_Cluster.pocket(i)

    def _gen_d_pockets(self) -> object:
        """
        Generate d-pocket dictionary of list of indices
        :return: dict of d_pockets
        """
        assert self._data is not None

        pockets = []
        lining_atoms = []
        for i in range(self.n_frames):
            for pocket in self.pockets(i):
                pockets.append(pocket)
                lining_atoms.append(pocket.lining_atoms)

        lining_atom_diff_matrix = np.zeros((len(lining_atoms),len(lining_atoms)))

        for i in range(len(lining_atoms)):
            for j in range(i,len(lining_atoms)):
                lining_atom_diff_matrix[i,j] = lining_atom_diff_matrix[j,i] = len(
                        lining_atoms[i].symmetric_difference(lining_atoms[j])) / len(
                    lining_atoms[i].union(lining_atoms[j]))

        clustered_list = list(fcluster(Z=linkage(squareform(lining_atom_diff_matrix),
                                                 method='complete'),
                                       t=self.config.dpocket_cluster_cutoff,
                                       criterion='distance'))

        d_pockets_idx_dict = {}
        for pocket,d_pocket_idx in zip(pockets,clustered_list):
            snapshot_idx = pocket.cluster.alpha_snapshot_idx
            pocket_idx = pocket.index
            if d_pocket_idx not in d_pockets_idx_dict:
                d_pockets_idx_dict[d_pocket_idx] = []
            d_pockets_idx_dict[d_pocket_idx].append((snapshot_idx,pocket_idx))

        for idx,item in d_pockets_idx_dict.items():
            yield AS_D_Pocket(self,pocket_indices=item)
