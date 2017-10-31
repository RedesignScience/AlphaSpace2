import numpy as np
from scipy.spatial import Voronoi, Delaunay
from scipy.cluster.hierarchy import linkage, fcluster
from mdtraj.core.topology import Topology, Residue, Atom
from mdtraj.core.trajectory import Trajectory
from mdtraj.core.element import *
from mdtraj import shrake_rupley
from AS_Funct import _tessellation, update_residue_method, update_atom_methods, getTetrahedronVolume, getContactMatrix


# noinspection PyAttributeOutsideInit,PyAttributeOutsideInit,PyAttributeOutsideInit,PyAttributeOutsideInit,PyAttributeOutsideInit,PyAttributeOutsideInit,PyTypeChecker
class AS_Cluster(object):
    def __init__(self, receptor, snapshot_idx=0):
        """
        Container for alpha, beta atoms and pocket
        :param receptor: AS_Struct
        :param snapshot_idx: int
        """
        self.snapshot_idx = snapshot_idx

        # Parent receptor attribute reference
        self.parent = receptor
        self.receptor_top = receptor.top

        self.universe = receptor.parent
        self.config = receptor.config

        # Initialize Container of topology and coordinate for pockets and alpha atoms
        self.top = Topology()

        self.top.add_chain()
        self.top.add_residue(name='ASC', chain=self.top.chain(0))
        self.traj = Trajectory(np.zeros([1, 0, 3]), self.top)

        # container for coordinate in located in the trajectory object, 0 frame
        self.xyz = self.traj.xyz

        # Link pocket name change to trajectory object

        # Initialize storage array
        self._contact = np.empty(self.n_alphas, int)
        self._polar_score = np.empty(self.n_alphas, float)
        self._nonpolar_score = np.empty(self.n_alphas, float)
        self._total_score = np.empty(self.n_alphas, float)

        # self._tessellation(self.config)

    def __call__(self, *args, **kwargs):
        self._tessellation(self.config)
        return self

    def __repr__(self):
        return "Alpha Atom cluster of #{} frame, {} pockets, {} Alpha Atoms".format(self.snapshot_idx, self.n_pockets,
                                                                                    self.n_alphas)

    @property
    def is_polar(self) -> np.array:
        return self.parent.is_polar

    @property
    def receptor_snapshot(self) -> object:
        return self.parent.trajectory[self.snapshot_idx]

    @property
    def n_alphas(self) -> int:
        return self.top.n_atoms

    @property
    def n_pockets(self) -> int:
        return self.top.n_residues

    @property
    def alphas(self) -> object:
        for atom in self.top._atoms:
            yield atom

    @property
    def pockets(self) -> object:
        return self.top.residues

    def pocket(self, index: int) -> object:
        return self.top.residue(index)

    def alpha(self, i):
        return self.top.atom(i)

    def residue(self, i):
        return self.top.residue(i)

    def _tessellation(self, config):
        """
        perform tessellation in order to generate the cluster of alpha atoms.
        :param config: object
        """
        # Generate Raw Tessellation simplices
        raw_alpha_lining_idx = Delaunay(self.receptor_snapshot.xyz[0]).simplices
        # Take coordinates from xyz file
        raw_alpha_lining_xyz = np.take(self.receptor_snapshot.xyz[0], raw_alpha_lining_idx[:, 0].flatten(), axis=0)

        # generate alpha atom coordinates
        raw_alpha_xyz = Voronoi(self.receptor_snapshot.xyz[0]).vertices

        # Calculate alpha sphere radii
        raw_alpha_sphere_radii = np.linalg.norm(raw_alpha_lining_xyz - raw_alpha_xyz, axis=1)

        # Filter the data based on radii cutoff
        filtered_alpha_idx = np.where(np.logical_and(config.min_r / 10.0 <= raw_alpha_sphere_radii,
                                                     raw_alpha_sphere_radii <= config.max_r / 10.0))[0]

        self.alpha_lining = np.take(raw_alpha_lining_idx, filtered_alpha_idx, axis=0)

        filtered_alpha_xyz = np.take(raw_alpha_xyz, filtered_alpha_idx, axis=0)

        # cluster the remaining vertices to assign index of belonging pockets
        zmat = linkage(filtered_alpha_xyz, method='average')
        pocket_idx = fcluster(zmat, self.config.pocket_cluster_distance / 10, criterion='distance')  # /10 turn A to nm

        self.alpha_pocket_index = pocket_idx - 1  # because cluster index start from 1

        # print(max(pocket_idx),min(pocket_idx))

        # Reorganize into list of pockets
        # self.pocket_alpha_atoms = [[] for _ in range(max(cluster))]
        # for alpha_cluster_i, alpha_atom_idx in sorted(
        #         zip(self.alpha_pocket_index, range(len(self.alpha_pocket_index)))):
        #     self.pocket_alpha_atoms[alpha_cluster_i].append(alpha_atom_idx)

        # Generate Residue container for pockets and add in atoms as AAC
        for _ in range(max(pocket_idx) - 1):
            residue = self.top.add_residue(name='ASC', chain=self.top.chain(0))
            residue.cluster = self
        for i, pocket_index in enumerate(self.alpha_pocket_index):
            atom = self.top.add_atom('AAC', None, self.top.residue(pocket_index), pocket_index)
            atom.index = i

        update_atom_methods(Atom)
        update_residue_method(Residue)

        for pocket in self.pockets:
            pocket.cluster = self  # assign the parent cluster
            alpha_index = [alpha.index for alpha in pocket.atoms]
            pocket.lining_atom_idx = self._get_lining_atoms(alpha_index)  # Assign pocket lining atoms
            pocket.lining_residue_idx = self._get_lining_residues(alpha_index)
        # Load trajectories
        self.traj.xyz = np.expand_dims(filtered_alpha_xyz, axis=0)
        filtered_lining_xyz = np.take(self.receptor_snapshot.xyz[0], self.alpha_lining, axis=0)
        # calculate the polarity of alpha atoms
        self._total_score = np.array([getTetrahedronVolume(i) for i in filtered_lining_xyz]) * 1000
        # here the 1000 is to convert nm^3 to A^3
        pocket_sasa = np.take(self._get_SASA(), self.alpha_lining)

        polar_ratio = np.average(np.take(self.is_polar, self.alpha_lining), axis=1, weights=pocket_sasa)

        self._polar_score = self._total_score * polar_ratio

        self._nonpolar_score = self._total_score - self._polar_score

        self._contact = np.zeros(self.top.n_atoms, dtype=int)


    def _get_SASA(self):
        """
        Calculate the absolute solvent accessible surface area.
        First calculate the SASA of the receptor by itself, then subtract it with sasa with the AAC.
        AAC are set to resemble Carbon with a radii - 0.17
        :return: np.array. The difference.
        """

        joined_traj = self.parent.trajectory[self.snapshot_idx].stack(self.traj)
        parent_sasa = shrake_rupley(self.parent.trajectory[self.snapshot_idx])[0]

        # todo fix the issue where too close of two alpha atoms cause shrake_rupley to crush

        joined_sasa = shrake_rupley(joined_traj, change_radii={'VS': 0.17})[0][:len(parent_sasa)]
        return parent_sasa - joined_sasa

    def _get_contact_list(self, binder_traj: object = None) -> np.array:
        """
        get list of alpha atom contact as bool
        :param binder_traj: object
        :return: np.array
        """
        if binder_traj is None:
            binder_traj = self.universe.binder.traj
        contact_matrix = getContactMatrix(self.traj.xyz[0], binder_traj.xyz[0],
                                          threshold=self.config.contact_threshold / 10)
        self._contact = np.array(np.any(contact_matrix, axis=1))
        return self._contact

    def _get_contact_space(self):
        """
        :return array
        """
        self._get_contact_list(self.universe.binder.traj)
        contact_space = self._contact * self._total_score
        return contact_space

    def _slice(self, new_alpha_indices: list) -> None:
        """
        This updates the alpha atom indexing after removing ones no in indices
        :param new_alpha_indices: list or set
        :return: None
        """
        if type(new_alpha_indices) != list:
            new_alpha_indices = list(new_alpha_indices)
        new_alpha_indices.sort()
        self.traj.atom_slice(new_alpha_indices, inplace=True)
        self.top = self.traj.top

        self.alpha_lining = np.take(self.alpha_lining, new_alpha_indices, axis=0)
        self._contact = np.take(self._contact, new_alpha_indices, axis=0)
        self._total_score = np.take(self._total_score, new_alpha_indices, axis=0)
        self._polar_score = np.take(self._polar_score, new_alpha_indices, axis=0)
        self._nonpolar_score = np.take(self._nonpolar_score, new_alpha_indices, axis=0)

        for pocket_idx, pocket in enumerate(self.pockets):
            pocket.cluster = self
            pocket.index = pocket_idx
        for atom_idx, atom in enumerate(self.top._atoms):
            atom.index = atom_idx

    def _get_lining_atoms(self, index_list: list) -> np.array:
        """
        get a set of surface lining atoms in the given cluster of alpha atoms
        :type index_list: alpha atom list
        """
        idx = np.take(self.alpha_lining, index_list, axis=0)
        total_list = idx.flatten()

        return np.unique(total_list)

    def _get_pocket_lining_atoms(self, pocket):
        return set(self._get_lining_atoms([int(atom.index) for atom in pocket.atoms]))

    def _get_lining_residues(self, index_list):
        # lining_atom_idx = np.take(self.alpha_atom_lining, index_list).flatten()

        lining_atom = [self.receptor_top.atom(i) for i in self._get_lining_atoms(index_list)]
        lining_residue_idx = [atom.residue.index for atom in lining_atom]

        return np.unique(lining_residue_idx)

    def _get_pocket_lining_residue(self, pocket):
        return self._get_lining_residues([atom.index for atom in pocket.alpha_atoms])

    def _get_cluster_centroid(self, cluster):
        xyz = np.take(self.traj.xyz[0], cluster, axis=0)
        return np.mean(xyz, axis=0)


class AS_D_Pocket:
    def __init__(self, parent_universe, pocket_indices=None):
        """
        Initialize a mask for information storage of a dpocket, this is user accessible and can only be read from public
         methods
        :param parent_universe: AS_Universe, parent universe
        """
        self.universe = parent_universe
        if pocket_indices is not None:
            self._pocket_indices = pocket_indices  # list of tuple, (snapshot_idx, pocket_idx)
        self.index = -1  # index of d-pocket in all d-pocket
        self._active = True

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
        Get pocket by index in d pocket
        :param key: int
        :return:
        """
        snapshot_idx, pocket_idx = self._pocket_indices[key]
        return self.universe.cluster(snapshot_idx).pocket(pocket_idx)

    def _index_to_pocket(self, index: tuple) -> object:
        return self.universe.cluster(snapshot_idx=index[0]).pocket(index[1])

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
        :param pocket: object, AS_Pocket
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
