import numpy as np
from scipy.spatial import Voronoi, Delaunay
from scipy.cluster.hierarchy import linkage, fcluster
from AS.AS_Funct import *
from mdtraj.core.topology import Topology, Residue, Atom
from mdtraj.core.trajectory import Trajectory
from mdtraj.core.element import *
from mdtraj import shrake_rupley


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
        self.is_polar = receptor.is_polar

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

        self._tessellation(self.config)

    def __repr__(self):
        return "Alpha Atom cluster of #{} frame, {} pockets, {} Alpha Atoms".format(self.snapshot_idx,self.n_pockets,self.n_alphas)

    @property
    def receptor_snapshot(self):
        return self.parent.trajectory[self.snapshot_idx]

    @property
    def n_alphas(self):
        return self.top.n_atoms

    @property
    def n_pockets(self):
        return self.top.n_residues

    @property
    def alphas(self):
        for atom in self.top._atoms:
            yield atom

    @property
    def pockets(self):
        return self.top.residues

    def pocket(self,index):
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
        for _ in range(max(pocket_idx)-1):
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
            pocket.lining_atom_idx = self._get_lining_atoms(alpha_index)   # Assign pocket lining atoms
            pocket.lining_residue_idx = self._get_lining_residues(alpha_index)
        # Load trajectories
        self.traj.xyz = np.expand_dims(filtered_alpha_xyz, axis=0)
        filtered_lining_xyz = np.take(self.receptor_snapshot.xyz[0], self.alpha_lining, axis=0)
        # calculate the polarity of alpha atoms
        self._total_score = np.array([getTetrahedronVolume(i) for i in filtered_lining_xyz])

        pocket_sasa = np.take(self._get_SASA(), self.alpha_lining)

        polar_ratio = np.average(np.take(self.is_polar, self.alpha_lining), axis=1, weights=pocket_sasa)

        self._polar_score = self._total_score * polar_ratio

        self._nonpolar_score = self._total_score - self._polar_score

        self._contact = np.zeros(self.top.n_atoms, dtype=int)

        # centoids = []
        # for pocket in self.pockets:
        #     centoid = np.zeros([3])
        #     for index in [alpha.index for alpha in pocket.atoms]:
        #         centoid = centoid + self.traj.xyz[0][index]
        #     centoids.append(centoid/pocket.n_atoms)
        #
        # n
        #
        # exit()


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

    def _get_contact_list(self, binder_traj=None):
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


        self.alpha_lining = np.take(self.alpha_lining, new_alpha_indices,axis=0)
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
        idx = np.take(self.alpha_lining,index_list,axis=0)
        total_list = idx.flatten()

        return np.unique(total_list)

    def _get_pocket_lining_atoms(self, pocket):
        return self._get_lining_atoms([atom.index for atom in pocket.atoms])

    def _get_lining_residues(self, index_list):
        # lining_atom_idx = np.take(self.alpha_atom_lining, index_list).flatten()

        lining_atom = [self.receptor_top.atom(i) for i in self._get_lining_atoms(index_list)]
        lining_residue_idx = [atom.residue.index for atom  in lining_atom]

        return np.unique(lining_residue_idx)

    def _get_pocket_lining_residue(self, pocket):
        return self._get_lining_residues([atom.index for atom in pocket.alpha_atoms])

    def _get_cluster_centroid(self, cluster):
        xyz = np.take(self.traj.xyz[0], cluster, axis=0)
        return np.mean(xyz, axis=0)


class AS_DPocket:
    def __init__(self, universe):
        self.universe = universe
        self._pocket_list = []

        self.index = -1

    def __getitem__(self, key):
        return self._pocket_list[key]

    def __contains__(self, pocket):
        return any(p == pocket for p in self._pocket_list)

    def __iter__(self):
        for i in self._pocket_list:
            yield i

    def add(self, pocket):
        """
        add a pocket to this dpocket
        :param pocket: object, AS_Pocket
        :return: bool, if the pocket isn't already in the dpocket list
        """
        if any(p == pocket for p in self._pocket_list):
            return False
        else:
            self._pocket_list.append(pocket)
            return True

    def pop(self, key):
        """
        pop a pocket
        :param key: int
        :return: AS_Pocket
        """
        return self._pocket_list.pop(key)

    def n_pockets(self):
        """
        total number of pockets
        :return: int
        """
        return len(self._pocket_list)

    def get_snapshots(self):
        """
        covered snapshot indices in a set
        :return: set
        """
        return iter(set([pocket.snapshot_idx for pocket in self]))

    def n_snapshot(self):
        """
        total number of covered snapshot
        :return:
        """
        return len(list(self.get_snapshots()))

    def pocket_center_xyz(self):
        """
        get the centroid of the each pockets
        :return: numpy array n_pockets * 3
        """
        return np.array([pocket.get_centroid for pocket in self])

    def alpha_atom_xyz(self):
        """
        get the xyz of each alpha_atoms in each pockets
        :return: np.array N * 3
        """
        return np.concatenate([pocket.xyz for pocket in self], axis=0)
