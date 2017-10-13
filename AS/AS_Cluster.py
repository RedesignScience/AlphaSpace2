import numpy as np
from scipy.spatial import Voronoi, Delaunay
from scipy.cluster.hierarchy import linkage, fcluster
from AS.AS_Funct import *
from mdtraj.core.topology import Topology
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
        self.receptor_snapshot = self.parent.trajectory[snapshot_idx]
        self.config = receptor.config
        self.is_polar = receptor.is_polar

        # Initialize Container of topology and coordinate for pockets and alpha atoms
        self.top = self.pockets = Topology()
        self.top.add_chain()
        self.top.add_residue(name='ASC', chain=self.top.chain(0))
        self.traj = Trajectory(np.zeros([1, 0, 3]), self.top)

        # container for coordinate in located in the trajectory object, 0 frame
        self.xyz = self.traj.xyz

        # Link pocket name change to trajectory object
        self.n_alpha = self.top.n_atoms
        self.n_pockets = self.top.n_residues
        self.alphas = self.top.atoms
        self.pockets = self.top.residues
        self.alpha = self.top.atom
        self.pocket = self.top.residue

        # Initialize storage array
        self.contact = np.empty(self.n_alpha, int)
        self.polar_score = np.empty(self.n_alpha, float)
        self.nonpolar_score = np.empty(self.n_alpha, float)
        self.total_score = np.empty(self.n_alpha, float)

        self.tessellation(self.config)
        if self.parent.structure_type == 0:
            self.oppsite_struct = self.parent.parent.binder
        elif self.parent.structure_type == 1:
            self.oppsite_struct = self.parent.parent.receptor
        else:
            self.oppsite_struct = None

    def __repr__(self):
        return "Alpha Atom cluster of #{} frame".format(self.snapshot_idx)

    def tessellation(self, config):
        """
        perform tessellation in order to generate the cluster of alpha atoms.
        :param config: object
        """
        # Generate Raw Tessellation vertices
        raw_alpha_lining = Delaunay(self.receptor_snapshot.xyz[0]).simplices
        raw_alpha_xyz = Voronoi(self.receptor_snapshot.xyz[0]).vertices
        raw_lining_xyz = np.take(self.receptor_snapshot.xyz[0], raw_alpha_lining[:, 0].flatten(), axis=0)

        # Calculate Raw alpha sphere radii
        raw_alpha_sphere_radii = np.linalg.norm(raw_lining_xyz - raw_alpha_xyz, axis=1)

        # Filter the data based on radii
        filtered_idx = np.where(np.logical_and(config.min_r / 10.0 <= raw_alpha_sphere_radii,
                                               raw_alpha_sphere_radii <= config.max_r / 10.0))[0]
        filtered_lining = np.take(raw_alpha_lining, filtered_idx, axis=0)

        self.alpha_atom_lining = filtered_lining

        filtered_xyz = np.take(raw_alpha_xyz, filtered_idx, axis=0)

        # cluster the remaining vertices to assign index of belonging pockets
        zmat = linkage(filtered_xyz, method='average')
        cluster = fcluster(zmat, self.config.pocket_cluster_distance / 10, criterion='distance')  # /10 turn A to nm

        self.alpha_atom_pocket_index = cluster - 1  # because cluster index start from 1

        # Reorganize into list of pockets
        self.pocket_alpha_atoms = [[] for _ in range(max(self.alpha_atom_pocket_index) + 1)]
        for alpha_cluster_i, alpha_atom_idx in sorted(
                zip(self.alpha_atom_pocket_index, range(len(self.alpha_atom_pocket_index)))):
            self.pocket_alpha_atoms[alpha_cluster_i].append(alpha_atom_idx)

        # Generate Residue container for pockets and add in atoms as AAC
        for _ in self.pocket_alpha_atoms[1:]:
            self.top.add_residue(name='ASC', chain=self.top.chain(0))
        for i, pocket_index in enumerate(self.alpha_atom_pocket_index):
            self.top.add_atom('AAC', None, self.top.residue(self.alpha_atom_pocket_index[i]), i)
        update_residue_method(self.pocket(0))
        update_atom_methods(self.alpha(0))

        for pocket in self.pockets:
            alpha_index = [alpha.index for alpha in pocket.alpha_atoms]
            pocket.lining_atom_idx = self._get_lining_atoms(alpha_index)
            pocket.lining_residue_idx = self._get_lining_residues(alpha_index)
        # Load trajectories
        self.traj.xyz = np.expand_dims(filtered_xyz, axis=0)
        filtered_lining_xyz = np.take(self.receptor_snapshot.xyz[0], self.alpha_atom_lining, axis=0)
        # calculate the polarity of alpha atoms
        self.total_score = np.array([getTetrahedronVolume(i) for i in filtered_lining_xyz])

        pocket_sasa = np.take(self._get_SASA(), self.alpha_atom_lining)

        polar_ratio = np.average(np.take(self.is_polar, self.alpha_atom_lining), axis=1, weights=pocket_sasa)

        self.polar_score = self.total_score * polar_ratio

        self.nonpolar_score = self.total_score - self.polar_score

        self.contact = np.zeros(self.top.n_atoms, dtype=int)

    def _get_SASA(self):
        """
        Calculate the absolute solvent accessible surface area.
        First calculate the sasa of the receptor by itself, then subtract it with sasa with the AAC.
        AAC are set to resemble Carbon with a radii - 0.17
        :return: np.array. The difference.
        """

        joined_traj = self.parent.trajectory[self.snapshot_idx].stack(self.traj)
        parent_sasa = shrake_rupley(self.parent.trajectory[self.snapshot_idx])[0]

        joined_sasa = shrake_rupley(joined_traj, change_radii={'VS': 0.17})[0][:len(parent_sasa)]
        return parent_sasa - joined_sasa

    def _get_contact_list(self, binder_traj):
        """

        :param binder_traj:
        :return:
        """
        self.contact = getContactMatrix(self.traj.xyz[0], binder_traj.xyz[0],
                                        threshold=self.config.contact_threshold / 10)
        return self.contact

    #
    # def _screen(self,index):
    #
    #     contact_alpha_atoms = np.where(self.alpha_atom_contact > 0)[0]
    #     contact_pocket = np.unique(np.take(self.alpha_atom_pocket_index, contact_alpha_atoms))
    #
    #     contact_pocket_alpha_atom_index = np.array(np.where(self.alpha_atom_pocket_index == contact_pocket[0]))
    #     for i in contact_pocket[1:]:
    #         contact_pocket_alpha_atom_index = np.append(contact_pocket_alpha_atom_index,
    #                                                     np.array(np.where(self.alpha_atom_pocket_index == i)))
    #     self.slice(contact_pocket_alpha_atom_index)

    def _get_contact_space(self):
        """
        calculate which point in the alpha cluster is in contact with the given coordinates
        :return array
        """
        self._get_contact_list(self.oppsite_struct.traj)
        return np.transpose(np.transpose(self.total_score) * self.contact)

    def _slice(self, index_list):
        """
        This updates the alpha atom indexing after removing noncontact ones
        :param index_list: list
        :return: None
        """
        self.traj.atom_slice(index_list, inplace=True)
        self.top = self.traj.top

        for pocket_idx, pocket in enumerate(self.pockets):
            pocket.index = pocket_idx
        for atom_idx, atom in enumerate(self.alphas):
            atom.index = atom_idx
        # self.alpha_atom_pocket_index = [atom.residue.index for atom in self.top.atoms]
        # self.alpha_atom_lining = np.take(self.alpha_atom_lining, index_list, axis=0)

        # update array storage
        self.contact = np.take(self.contact, index_list, axis=0)
        self.total_score = np.take(self.total_score, index_list, axis=0)
        self.polar_score = np.take(self.polar_score, index_list, axis=0)
        self.nonpolar_score = np.take(self.nonpolar_score, index_list, axis=0)

    def _get_lining_atoms(self, index_list):
        """
        get a set of surface lining atoms in the given cluster of alpha atoms
        :type index_list: list
        """
        return np.unique(np.take(self.alpha_atom_lining, index_list).flatten())

    def _get_pocket_lining_atoms(self, pocket):

        return self._get_lining_atoms([atom.index for atom in pocket.alpha_atoms])

    def _get_lining_residues(self, index_list):
        # lining_atom_idx = np.take(self.alpha_atom_lining, index_list).flatten()

        lining_atom = [self.receptor_top.atom(i) for i in self._get_lining_atoms(index_list)]
        lining_residue_idx = [atom.residue.index for atom in lining_atom]

        return np.unique(lining_residue_idx)

    def _get_pocket_lining_residue(self, pocket):
        return self._get_lining_residues([atom.index for atom in pocket.alpha_atoms])

    def get_cluster_centroid(self, cluster):
        xyz = np.take(self.traj.xyz[0], cluster, axis=0)
        return np.mean(xyz, axis=0)

#
# class AS_Alpha:
#     def __init__(self, atom, cluster, contact=False):
#         self.atom_object = atom
#         self.contact = contact
#         self.cluster = cluster
#         self.index = self.atom_object.index
#
#     def get_index(self):
#         return self.atom_object.index
#
#     def get_pocket(self):
#         return AS_Pocket(self.atom_object.residue, cluster=self.cluster)
#
#     def get_polar_score(self):
#         return self.cluster.polar_score[self.index]
#
#     def get_nonpolar_score(self):
#         return self.cluster.nonpolar_score[self.index]
#
#     def get_contact(self):
#         return self.contact


# class AS_Pocket:
#     def __init__(self, contact_score, total_score, parent_cluster, index=0):
#         self.index = index
#         self.parent_cluster = parent_cluster
#         self.snapshot_idx = self.parent_cluster.snapshot_idx
#         self._contact_score = contact_score
#         self._total_score = total_score
#         self._component_lookup = {'all': 0, 'polar': 1, 'nonpolar': 2}
#
#         self.alpha_atoms = self.parent_cluster.pocket[self.index].atoms
#         self.parent_dpocket = None
#
#         self.lining_atom_index = set()
#
#         self.is_interface = False
#
#     def __eq__(self, as_pocket):
#         return (self.snapshot_idx, self.index) == (self.snapshot_idx, self.index)
#
#     def xyz(self):
#         return self.parent_cluster.xyz[[atom.index for atom in self.alpha_atoms]]
#
#     def get_centroid(self):
#         return self.parent_cluster.get_cluster_centroid([atom.index for atom in self.alpha_atoms])
#
#     def get_contact_score(self, component='all'):
#         if component in self._component_lookup:
#             return self._contact_score[self._component_lookup[component]]
#         else:
#             raise ("Invalid component name, please use all, polar, or nonpolar, you used {}".format(component))
#
#     def get_total_score(self, component='all'):
#         if component in self._component_lookup:
#             return self._total_score[self._component_lookup[component]]
#         else:
#             raise ("Invalid component name, please use all, polar, or nonpolar, you used {}".format(component))
#
#     def __repr__(self):
#         return "Pocket {} with total score of {}".format(self.index, self.get_total_score('all'))
#
#     def get_lining_atom_index(self):
#         if not self.lining_atom_index:
#             self.lining_atom_index = self.parent_cluster._get_pocket_lining_atoms(self)
#
#         return self.lining_atom_index


class AS_DPocket:
    def __init__(self, parent_session):
        self.session = parent_session
        self.pocket_list = []
        self.index = -1

    def __getitem__(self, item):
        return self.pocket_list[item]

    def __contains__(self, pocket):
        return any(p == pocket for p in self.pocket_list)

    def __iter__(self):
        for i in self.pocket_list:
            yield i

    def add(self, pocket):
        """
        add a pocket to this dpocket
        :param pocket: object, AS_Pocket
        :return: bool, if the pocket isn't already in the dpocket list
        """
        if any(p == pocket for p in self.pocket_list):
            return False
        else:
            self.pocket_list.append(pocket)
            return True

    def pop(self, key):
        """
        pop a pocket
        :param key: int
        :return: AS_Pocket
        """
        return self.pocket_list.pop(key)

    def n_pockets(self):
        """
        total number of pockets
        :return: int
        """
        return len(self.pocket_list)

    def get_snapshots(self):
        """
        covered snapshot indices
        :return: set
        """
        return set([pocket.snapshot_idx for pocket in self])

    def n_snapshot(self):
        """
        total number of covered snapshot
        :return:
        """
        return len(self.get_snapshots())

    def pocket_xyz(self):
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
