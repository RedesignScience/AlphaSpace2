import numpy as np


# noinspection PyAttributeOutsideInit,PyAttributeOutsideInit,PyAttributeOutsideInit,PyAttributeOutsideInit,PyAttributeOutsideInit,PyAttributeOutsideInit,PyTypeChecker
# class AS_Cluster(object):
#     def __init__(self, receptor, snapshot_idx=0):
#         """
#         Container for alpha, beta atoms and pocket
#         :param receptor: AS_Struct
#         :param snapshot_idx: int
#         """
#         self.snapshot_idx = snapshot_idx
#
#         # Parent receptor attribute reference
#         self.parent = receptor
#         self.receptor_top = receptor.top
#
#         self.universe = receptor.parent
#         self.config = receptor.config
#
#         # Initialize Container of topology and coordinate for pockets and alpha atoms
#         self.top = Topology()
#
#         self.top.add_chain()
#         self.top.add_residue(name='ASC', chain=self.top.chain(0))
#         self.traj = Trajectory(np.zeros([1, 0, 3]), self.top)
#
#         # container for coordinate in located in the trajectory object, 0 frame
#         self.xyz = self.traj.xyz
#
#         # Link pocket name change to trajectory object
#
#         # Initialize storage array
#         self._contact = np.empty(self.n_alphas, int)
#         self._polar_score = np.empty(self.n_alphas, float)
#         self._nonpolar_score = np.empty(self.n_alphas, float)
#         self._total_score = np.empty(self.n_alphas, float)
#
#         # self._tessellation(self.config)
#
#     def __call__(self, *args, **kwargs):
#         self._tessellation(self.config)
#         return self
#
#     def __repr__(self):
#         return "Alpha Atom cluster of #{} frame, {} pockets, {} Alpha Atoms".format(self.snapshot_idx, self.n_pockets,
#                                                                                     self.n_alphas)
#
#
#     def n_alphas(self) -> int:
#         return self.top.n_atoms
#
#     @property
#     def n_pockets(self) -> int:
#         return self.top.n_residues
#
#     @property
#     def alphas(self) -> object:
#         for atom in self.top._atoms:
#             yield atom
#
#     @property
#     def pockets(self) -> object:
#         return self.top.residues
#
#     def pocket(self, index: int) -> object:
#         return self.top.residue(index)
#
#     def alpha(self, i):
#         return self.top.atom(i)
#
#     def residue(self, i):
#         return self.top.residue(i)
#
#     def _tessellation(self, config):
#         """
#         perform tessellation in order to generate the cluster of alpha atoms.
#         :param config: object
#         """
#         # Generate Raw Tessellation simplices
#         raw_alpha_lining_idx = Delaunay(self.receptor_snapshot.xyz[0]).simplices
#         # Take coordinates from xyz file
#         raw_alpha_lining_xyz = np.take(self.receptor_snapshot.xyz[0], raw_alpha_lining_idx[:, 0].flatten(), axis=0)
#
#         # generate alpha atom coordinates
#         raw_alpha_xyz = Voronoi(self.receptor_snapshot.xyz[0]).vertices
#
#         # Calculate alpha sphere radii
#         raw_alpha_sphere_radii = np.linalg.norm(raw_alpha_lining_xyz - raw_alpha_xyz, axis=1)
#
#         # Filter the data based on radii cutoff
#         filtered_alpha_idx = np.where(np.logical_and(config.min_r / 10.0 <= raw_alpha_sphere_radii,
#                                                      raw_alpha_sphere_radii <= config.max_r / 10.0))[0]
#
#         self.alpha_lining = np.take(raw_alpha_lining_idx, filtered_alpha_idx, axis=0)
#
#         filtered_alpha_xyz = np.take(raw_alpha_xyz, filtered_alpha_idx, axis=0)
#
#         # cluster the remaining vertices to assign index of belonging pockets
#         zmat = linkage(filtered_alpha_xyz, method='average')
#         pocket_idx = fcluster(zmat, self.config.pocket_cluster_distance / 10, criterion='distance')  # /10 turn A to nm
#
#         self.alpha_pocket_index = pocket_idx - 1  # because cluster index start from 1
#
#         # print(max(pocket_idx),min(pocket_idx))
#
#         # Reorganize into list of pockets
#         # self.pocket_alpha_atoms = [[] for _ in range(max(cluster))]
#         # for alpha_cluster_i, alpha_atom_idx in sorted(
#         #         zip(self.alpha_pocket_index, range(len(self.alpha_pocket_index)))):
#         #     self.pocket_alpha_atoms[alpha_cluster_i].append(alpha_atom_idx)
#
#         # Generate Residue container for pockets and add in atoms as AAC
#         for _ in range(max(pocket_idx) - 1):
#             residue = self.top.add_residue(name='ASC', chain=self.top.chain(0))
#             residue.cluster = self
#         for i, pocket_index in enumerate(self.alpha_pocket_index):
#             atom = self.top.add_atom('AAC', None, self.top.residue(pocket_index), pocket_index)
#             atom.index = i
#
#         update_atom_methods(Atom)
#         update_residue_method(Residue)
#
#         for pocket in self.pockets:
#             pocket.cluster = self  # assign the parent cluster
#             alpha_index = [alpha.index for alpha in pocket.atoms]
#             pocket.lining_atom_idx = self._get_lining_atoms(alpha_index)  # Assign pocket lining atoms
#             pocket.lining_residue_idx = self._get_lining_residues(alpha_index)
#         # Load trajectories
#         self.traj.xyz = np.expand_dims(filtered_alpha_xyz, axis=0)
#         filtered_lining_xyz = np.take(self.receptor_snapshot.xyz[0], self.alpha_lining, axis=0)
#         # calculate the polarity of alpha atoms
#         self._total_score = np.array([getTetrahedronVolume(i) for i in filtered_lining_xyz]) * 1000
#         # here the 1000 is to convert nm^3 to A^3
#         pocket_sasa = np.take(self._get_SASA(), self.alpha_lining)
#
#         polar_ratio = np.average(np.take(self.is_polar, self.alpha_lining), axis=1, weights=pocket_sasa)
#
#         self._polar_score = self._total_score * polar_ratio
#
#         self._nonpolar_score = self._total_score - self._polar_score
#
#         self._contact = np.zeros(self.top.n_atoms, dtype=int)
#
#     def _get_SASA(self):
#         """
#         Calculate the absolute solvent accessible surface area.
#         First calculate the SASA of the receptor by itself, then subtract it with sasa with the AAC.
#         AAC are set to resemble Carbon with a radii - 0.17
#         :return: np.array. The difference.
#         """
#
#         joined_traj = self.parent.trajectory[self.snapshot_idx].stack(self.traj)
#         parent_sasa = shrake_rupley(self.parent.trajectory[self.snapshot_idx])[0]
#
#         # todo fix the issue where too close of two alpha atoms cause shrake_rupley to crush
#
#         joined_sasa = shrake_rupley(joined_traj, change_radii={'VS': 0.17})[0][:len(parent_sasa)]
#         return parent_sasa - joined_sasa
#
#     def _get_contact_list(self, binder_traj: object = None) -> np.array:
#         """
#         get list of alpha atom contact as bool
#         :param binder_traj: object
#         :return: np.array
#         """
#         if binder_traj is None:
#             binder_traj = self.universe.binder.traj
#         contact_matrix = getContactMatrix(self.traj.xyz[0], binder_traj.xyz[0],
#                                           threshold=self.config.contact_threshold / 10)
#         self._contact = np.array(np.any(contact_matrix, axis=1))
#         return self._contact
#
#     def _get_contact_space(self):
#         """
#         :return array
#         """
#         self._get_contact_list(self.universe.binder.traj)
#         contact_space = self._contact * self._total_score
#         return contact_space
#
#     def _slice(self, new_alpha_indices: list) -> None:
#         """
#         This updates the alpha atom indexing after removing ones no in indices
#         :param new_alpha_indices: list or set
#         :return: None
#         """
#         if type(new_alpha_indices) != list:
#             new_alpha_indices = list(new_alpha_indices)
#         new_alpha_indices.sort()
#         self.traj.atom_slice(new_alpha_indices, inplace=True)
#         self.top = self.traj.top
#
#         self.alpha_lining = np.take(self.alpha_lining, new_alpha_indices, axis=0)
#         self._contact = np.take(self._contact, new_alpha_indices, axis=0)
#         self._total_score = np.take(self._total_score, new_alpha_indices, axis=0)
#         self._polar_score = np.take(self._polar_score, new_alpha_indices, axis=0)
#         self._nonpolar_score = np.take(self._nonpolar_score, new_alpha_indices, axis=0)
#
#         for pocket_idx, pocket in enumerate(self.pockets):
#             pocket.cluster = self
#             pocket.index = pocket_idx
#         for atom_idx, atom in enumerate(self.top._atoms):
#             atom.index = atom_idx
#
#     def _get_lining_atoms(self, index_list: list) -> np.array:
#         """
#         get a set of surface lining atoms in the given cluster of alpha atoms
#         :type index_list: alpha atom list
#         """
#         idx = np.take(self.alpha_lining, index_list, axis=0)
#         total_list = idx.flatten()
#
#         return np.unique(total_list)
#
#     def _get_pocket_lining_atoms(self, pocket):
#         return set(self._get_lining_atoms([int(atom.index) for atom in pocket.atoms]))
#
#     def _get_lining_residues(self, index_list):
#         # lining_atom_idx = np.take(self.alpha_atom_lining, index_list).flatten()
#
#         lining_atom = [self.receptor_top.atom(i) for i in self._get_lining_atoms(index_list)]
#         lining_residue_idx = [atom.residue.index for atom in lining_atom]
#
#         return np.unique(lining_residue_idx)
#
#     def _get_pocket_lining_residue(self, pocket):
#         return self._get_lining_residues([atom.index for atom in pocket.alpha_atoms])
#
#     def _get_cluster_centroid(self, cluster):
#         xyz = np.take(self.traj.xyz[0], cluster, axis=0)
#         return np.mean(xyz, axis=0)


class AS_D_Pocket:
    def __init__(self, parent_structure, pocket_indices=None):
        """
        Initialize a mask for information storage of a dpocket, this is user accessible and can only be read from public
         methods
        :param parent_structure: AS_Universe, parent universe
        """
        self.structure = parent_structure
        if pocket_indices is not None:
            self._pocket_indices = pocket_indices  # list of tuple, (snapshot_idx, pocket_idx)
        self.index = -1  # index of d-pocket in all d-pocket
        self._active = True
        self._data = self.structure._data

    @property
    def is_active(self):
        return self._active

    def activate(self):
        self._active = True

    def deactivate(self):
        self._active = False

    def __bool__(self):
        return self._active

    @property
    def _pocket_set(self):
        return set(self._pocket_indices)

    def __getitem__(self, key):
        """
        :param key:
        :return:
        """

    def _index_to_pocket(self, index: tuple) -> object:
        return self.structure.cluster(snapshot_idx=index[0]).pocket(index[1])

    def _pocket_to_index(self, pocket: object) -> tuple:
        return pocket.cluster.snapshot_index, pocket.index

    def __contains__(self, item) -> bool:
        if type(item) is not tuple:
            item = self._pocket_to_index(item)
        return any(p == item for p in self._pocket_indices)

    @property
    def __iter__(self) -> object:
        for idx in self._pocket_indices:
            yield self._index_to_pocket(idx)

    @property
    def pocket_indices(self) -> tuple:
        for i in self._pocket_indices:
            yield i

    @property
    def pockets(self):
        for pocket_idx in self.pocket_indices:
            yield self._index_to_pocket(pocket_idx)

    @property
    def pockets_by_snapshots(self):
        snapshot_dict = {}
        for ss_idx, pkt_idx in self._pocket_indices:
            if ss_idx not in snapshot_dict:
                snapshot_dict[ss_idx] = []
            snapshot_dict[ss_idx].append(self._index_to_pocket((ss_idx, pkt_idx)))
        return snapshot_dict

    def add(self, item):
        """
        add a pocket to this dpocket
        :param item: object, AS_Pocket
        :return: bool, if the pocket isn't already in the dpocket list
        """
        if type(item) is not tuple:
            item = self._pocket_to_index(item)
        if item in self:
            raise Exception('Pocket {} already exist in Frame {}'.format(item[1], item[0]))

    def pop(self, key):
        """
        pop a pocket
        :param key: int
        :return: AS_Pocket
        """
        return self._index_to_pocket(self._pocket_indices.pop(key))

    @property
    def n_pockets(self) -> int:
        """
        total number of pockets
        :return: int
        """
        return len(self._pocket_indices)

    @property
    def n_snapshots(self) -> int:
        """
        total number of covered snapshot
        :return:
        """
        return len(self.snapshots)

    @property
    def snapshots(self) -> list:
        """
        Return snapshot index for all the pockets in this d-pocket
        :return: list of int
        """
        return sorted(list(set([i[1] for i in self._pocket_indices])))

    @property
    def pocket_scores(self) -> tuple:
        """
        :return: iter of tuple (polar, nonpolar)
        """
        for pocket in self.pockets:
            yield pocket.get_polar_score(), pocket.get_nonpolar_score()

    @property
    def pocket_xyz(self) -> np.array:
        for pocket in self.pockets:
            yield pocket.get_centroid()


class AS_Pocket:
    """
    This is the pocket container for the topological information of a pocket
    A pocket is a collection of alpha atoms clustered by average linkage
    """

    def __init__(self, alpha_idx, snapshot_idx, idx, parent_structure):
        self._idx = idx
        self._snapshot_idx = snapshot_idx
        self._alpha_idx = alpha_idx
        self.parent_structure = parent_structure

    def __len__(self) -> int:
        return len(self._alpha_idx)

    # def merge_pocket(self,other):
    #     """
    #     Merge two pockets, the new index of the pocket will be the same as the first one
    #     :param other: AS_Pocket
    #     :return: AS_Pocket
    #     """
    #     if self._snapshot_idx != other._snapshot_idx:
    #         raise Exception("Can't merge pockets from two snapshots")
    #     return AS_Pocket(set(self._alpha_idx).union(other._alpha_idx),self._snapshot_idx,self._idx,self.parent_structure)
    #
    # def add_alpha(self,idx):
    #     if idx in self._alpha_idx:
    #         return False
    #     else:
    #         self._alpha_idx.append(idx)
    #         return self
    #
    # def __add__(self, other):
    #     if type(other) == AS_Pocket:
    #         return self.merge_pocket(other)
    #     elif type(self) == AS_AlphaAtom:
    #         return self.add_alpha(other.idx)


    @property
    def alphas(self):
        for i in self.alpha_idx:
            yield AS_AlphaAtom(i, self.parent_structure, self)

    @property
    def alpha_idx(self):
        return np.array(self._alpha_idx, dtype=int)

    @property
    def _data(self):
        return self.parent_structure._data

    @property
    def is_active(self, anyone=True):
        if anyone:
            return self._data.is_active(self.alpha_idx).any()
        else:
            return self._data.is_active(self.alpha_idx).all()

    @property
    def xyz(self) -> np.ndarray:
        return self._data.xyz(self.alpha_idx)

    def polar_score(self, decimals=0) -> float:
        return float(np.around(np.sum(self._data.polar_score(self.alpha_idx)), decimals=decimals))

    def polar_scores(self, decimals=0) -> np.ndarray:
        return np.around(self._data.polar_score(self.alpha_idx), decimals=decimals)

    def nonpolar_score(self, decimals=0) -> float:
        return float(np.around(np.sum(self._data.nonpolar_score(self.alpha_idx)), decimals=decimals))

    def nonpolar_scores(self, decimals=0) -> np.ndarray:
        return np.around(self._data.nonpolar_score(self.alpha_idx), decimals=decimals)

    def total_score(self, decimals=0):
        return self.polar_score(decimals=decimals) + self.nonpolar_score(decimals=decimals)

    @property
    def lining_atoms_idx(self):
        return np.unique(self._data.lining_atoms_idx(self.alpha_idx))

    @property
    def lining_atoms(self):
        for idx in self.lining_atoms_idx:
            yield self.parent_structure.top.atom(idx)

    @property
    def lining_residues_idx(self):
        return np.unique([atom.residue.index for atom in self.lining_atoms])

    @property
    def lining_residues(self) -> np.ndarray:
        residue_idx = self.lining_residues_idx
        for idx in residue_idx:
            yield self.parent_structure.top.residue(idx)

    def activate(self):
        self._data.activate(self.alpha_idx)

    def deactivate(self):
        self._data.deactivate(self.alpha_idx)


class AS_AlphaAtom:
    """
    This is the alpha atom container, which is the constituents of the pocket object.
    This object only contains reference of alpha atom id, and all info is stored in the array you can access by methods.
    """

    def __init__(self, idx, parent_structure, parent_pocket=None):
        self._idx = idx
        self.parent_structure = parent_structure
        self.parent_pocket = parent_pocket

    def __repr__(self):
        return "AlphaAtom"

    def __bool__(self):
        return self.is_active

    def __sub__(self, other):
        return np.linalg.norm(self.xyz - other.xyz)

    @property
    def _data(self) -> object:
        return self.parent_structure._data

    @property
    def is_active(self):
        return

    @property
    def idx(self):
        return self._idx

    @property
    def snapshot_idx(self):
        return self._data.snapshot_idx(self._idx)

    @property
    def xyz(self):
        return self._data.xyz(self._idx)

    @property
    def polar_score(self):
        return self._data.polar_score(self._idx)

    @property
    def nonpolar_score(self):
        return self._data.nonpolar_score(self._idx)

    @property
    def lining_atoms_idx(self):
        return self._data.lining_atoms_idx(self._idx)

    @property
    def is_contact(self):
        return self._data.is_contact(self._idx)


class AS_Data(np.ndarray):
    """
    Container object inherited from numpy array object, you can access information directly here

    Column, Content,   default value
    0       idx
    1       snapshot_idx
    2       x
    3       y
    4       z
    5       lining_atom_idx_1
    6       lining_atom_idx_1
    7       lining_atom_idx_1
    8       lining_atom_idx_1
    9       polar_score        0
    10      nonpolar_score     0
    11      is_active          1
    12      is_contact         0
    13      pocket_idx         0
    14      radii
    """

    def __new__(cls, array_data, parent_structure=None):
        obj = super(AS_Data, cls).__new__(cls, buffer=array_data, shape=array_data.shape)
        obj.parent_structure = parent_structure
        return obj

    def idx(self, idx):
        return self[np.array(idx), 0]

    def alpha_snapshot_idx(self, idx):
        return self[np.array(idx), 1]

    @property
    def snapshots_idx(self):
        return np.sort(np.unique(self[:, 1])).astype(int)

    def snapshot(self, snapshot_idx):
        return self[self.snapshot_alpha_idx(snapshot_idx)]

    def xyz(self, idx):
        return self[np.array(idx), 2:5]

    def polar_score(self, idx):
        return self[np.array(idx), 9]

    def nonpolar_score(self, idx):
        return self[np.array(idx), 10]

    def lining_atoms_idx(self, idx):
        return self[np.array(idx, dtype=int), 5:9].astype(int)

    def is_active(self, idx):
        return self[np.array(idx), 11]

    def is_contact(self, idx):
        return self[np.array(idx), 12]

    def pocket(self, idx):
        return self[np.array(idx), 13].astype(int)

    def snapshot_pockets_idx(self, snapshot_idx=0):
        return self[self[:, 1] == snapshot_idx][:, 13].astype(int)

    def snapshot_alpha_idx(self, snapshot_idx=0):
        return self[self[:, 1] == snapshot_idx][:, 0].astype(int)

    def radii(self, idx):
        return self[np.array(idx), 14]

    def activate(self, idx):
        self[idx, 11] = 1

    def deactivate(self, idx):
        self[idx, 11] = 0

    @property
    def all_idx(self):
        return self[:, 0].astype(int)
