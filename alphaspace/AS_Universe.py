from itertools import chain

import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
import mdtraj

from .AS_Cluster import AS_D_Pocket
from .AS_Config import AS_Config
from .AS_Funct import getCosAngleBetween, combination_intersection_count, combination_union_count
from .AS_Struct import AS_Structure


# noinspection PyAttributeOutsideInit,PyAttributeOutsideInit,PyAttributeOutsideInit,PyAttributeOutsideInit,PyAttributeOutsideInit,PyAttributeOutsideInit
class AS_Universe(object):
    def __init__(self, receptor=None, binder=None, guess_receptor_binder=True, guess_by_order=True, config=None,
                 tag=""):
        """
        Container for an AlphaSpace session, have child container receptor and binder
        :param receptor: object
        :param binder: object
        :param guess_by_order: book
        :param guess_receptor_binder: bool, guess is based on molecule size
        :param config: object,AS_config
        """

        self.config = config if config else AS_Config()

        self.set_receptor(receptor)
        self.set_binder(binder)

        if guess_receptor_binder and receptor and not binder:
            if self.guess_receptor_binder(receptor, guess_by_order):
                pass
            else:
                raise Exception('No binder detected')

        self._pocket_list = {}

        self.others = None
        self._view = None

        self._d_pockets = {}
        self._pocket_network = {}

        self.tag = tag




    def __repr__(self):
        return "Receptor of {} residues {} atoms | Binder of {} residues {} atoms".format(self.receptor.n_residues,
                                                                                          self.receptor.n_atoms,
                                                                                          self.binder.n_residues,
                                                                                          self.binder.n_atoms)

    @property
    def _data(self):
        return self.receptor._data

    @property
    def frames(self):
        for i in range(self.n_frames):
            yield i

    @property
    def n_frames(self):
        if self.receptor is None:
            return 0
        else:
            return self.receptor.trajectory.n_frames

    @property
    def n_atoms(self) -> int:
        """
        return the total number of atoms in receptor and binders
        :return: int
        """
        return self.receptor.n_atoms + self.binder.n_atoms

    @property
    def n_residues(self):
        """
        return the total number of residues in the receptor and binder
        :return: int
        """
        return self.receptor.n_residues + self.binder.n_residues

    @property
    def molecules(self):
        """
        iterate over receptor and binder, if there is any
        :return: iter
        """
        for m in self.receptor, self.binder:
            if m is not None:
                yield m
            else:
                continue

    @property
    def d_pockets(self) -> AS_D_Pocket:
        """
        calculate the d pockets and give a iterator os AS_D_Pocket
        :return: object, AS_D_Pocket
        """
        if not self._d_pockets:
            self._d_pockets = list(self.receptor._gen_d_pockets())
        for p in self._d_pockets:
            yield p

    def d_pocket(self, i) -> AS_D_Pocket:
        if not self._d_pockets:
            self._d_pockets = self.receptor._gen_d_pockets()
        return self._d_pockets[i]

    # def _is_processed(self, snapshot_idx: int) -> bool:
    #     if snapshot_idx in self.receptor._clusters:
    #         return True
    #     else:
    #         return False

    def pockets(self, snapshot_idx: int = 0, active_only: bool = True) -> list:
        if snapshot_idx not in self._pocket_list:
            self._pocket_list[snapshot_idx] = self.pocket_list(snapshot_idx, active_only)

        return self._pocket_list[snapshot_idx]

    def alphas(self, snapshot_idx=0):
        for pocket in self.receptor.pockets(snapshot_idx):
            for alpha in pocket.alphas:
                yield alpha

    def pocket(self, pocket_idx, snapshot_idx=0):
        return self.receptor.pocket(pocket_idx, snapshot_idx)

    def pocket_list(self, snapshot_idx: int = 0, active_only: bool = True):
        pocket_list = []
        # i = 1
        for pocket in self.receptor.pockets(snapshot_idx):
            if pocket.is_active or (not active_only):
                # pocket._reordered_index = i
                pocket_list.append(pocket)
            else:
                continue
        pocket_list.sort(key=lambda p: p.get_space(), reverse=True)
        i = 1
        for pocket in pocket_list:
            pocket._reordered_index = i
            i += 1

        return pocket_list




    def guess_receptor_binder(self, traj, by_order: bool = True) -> bool:
        """
        Divide receptor trajectory based on connectivity, set larger molecule as receptor.
        This process automatically removes water and other solvents
        :param traj: target trajectory
        :param by_order: bool, if False, guess by appearances in file
        :return: bool, if any macro molecule were found.
        """
        if traj is None:
            raise Exception("Cannot guess receptor and binder, no structure detected")
        molecule_list = []
        # remove solvent and divide into molecules. This will guess which one is the
        for molecule in traj.topology.find_molecules():
            if len(molecule) > 1 and not next(iter(molecule)).residue.is_water:
                molecule_list.append(molecule)
                # print(len(molecule))
            else:
                # print(len(molecule))
                continue

        # print(len(molecule_list))

        if not by_order:
            molecule_list.sort(key=len, reverse=True)
        if len(molecule_list) > 1:
            self.set_receptor(traj.atom_slice(np.sort(np.array([atom.index for atom in molecule_list[0]], dtype=int))))
            self.set_binder(traj.atom_slice(np.sort(np.array([atom.index for atom in molecule_list[1]], dtype=int))))
            return True
        elif len(molecule_list) == 1:
            self.set_receptor(traj.atom_slice([atom.index for atom in molecule_list[0]]))
            self.binder = None
            return False
        else:
            return False

    def set_binder(self, structure, append=False):
        """
        set binder (ligand) in session
        :param structure: object, trajectory
        :param append: Bool, if the new binder should be appended to the preview one, default overwritten.
        :return:
        """
        if structure is None:
            self.binder = None
            return
        if append and (self.binder is not None):
            x = self.binder.trajectory + structure
            self.binder.trajectory = x
        else:
            self.binder = AS_Structure(structure, structure_type=1, parent=self)

    def set_receptor(self, structure, append=False, keepH=True):
        """
        set receptor (protein) in session
        :param structure: trajectory
        :param append: Bool, if the new binder should be appended to the preview one, default overwritten.
        :return:
        """

        if structure is None:
            self.receptor = None
            return

        if not keepH:
            non_h_idx = structure.topology.select_atom_indices(selection='heavy')
            structure.atom_slice(non_h_idx, inplace=True)

        if append and (self.receptor is not None):
            x = self.receptor.trajectory + structure
            self.receptor.trajectory = x
        else:
            self.receptor = AS_Structure(structure, structure_type=0, parent=self)

    def run_alphaspace(self):
        self.run_alphaspace_mp()

    def run_alphaspace_mp(self, cpu=None):

        """
        :type cpu: int
        """
        import multiprocessing as mp
        from .AS_Funct import _tessellation
        from concurrent.futures import ProcessPoolExecutor

        if not self.binder:
            self.config.screen_by_face = self.config.screen_by_lig_cntct = False

        # Disable screen by contact and face when there is no binder.



        cpu = mp.cpu_count() if not cpu else int(cpu)

        def executor(argslist):
            with ProcessPoolExecutor(max_workers=cpu) as ex:
                return ex.map(_tessellation, argslist)

        receptor_ss = [self.receptor.traj[i] for i in range(self.n_frames)]
        if self.binder:
            binder_ss = [self.binder.traj[i] for i in range(self.n_frames)]
        else:
            binder_ss = [[] for i in range(self.n_frames)]
        config = [self.config for i in range(self.n_frames)]
        snapshot_indices = range(self.n_frames)
        is_polar = [self.receptor.is_polar for _ in range(self.n_frames)]
        arglist = list(zip(receptor_ss, binder_ss, config, snapshot_indices, is_polar))

        data_list = list(executor(arglist))

        self.receptor._combine_data(data_list)


    def _get_face_atoms(self):

        """
        Calculate the snapshot interface atom.
        The interface atom is defined as whose ASA is reduced with introduction of ligand.
        :return: ndarray num_snapshots * num_atom_receptor
        """
        receptor_snapshot = self.receptor.traj

        complex_snapshot = receptor_snapshot.stack(self.binder.traj)

        receptor_snapshot_sasa = mdtraj.shrake_rupley(receptor_snapshot)

        complex_snapshot_sasa = mdtraj.shrake_rupley(complex_snapshot)

        sasa_diff = receptor_snapshot_sasa - complex_snapshot_sasa[:, :self.receptor.n_atoms]

        return sasa_diff > 0

        # return np.where((sasa_diff > 0).any(axis=0))[0]

    def _gen_communities_legacy(self):
        import networkx
        self.communities = {}
        for snapshot_idx in range(self.n_frames):
            pocket_graph = networkx.Graph()
            pocket_graph.add_nodes_from(self.pockets(snapshot_idx, active_only=False))

            # def connectPockets(p1, p2):
            #     pocket_vector1 = p1.lining_atoms_centroid - p1.centroid
            #     pocket_vector2 = p2.lining_atoms_centroid - p2.centroid
            #     if getCosAngleBetween(pocket_vector1, pocket_vector2) > 0:  # pocket vector facing inwards
            #         pocket_graph.add_edge(p1, p2)
            #         return True
            #     return False

            contact_pair = np.array(np.where(combination_intersection_count(
                [pocket.lining_atoms_idx for pocket in self.pockets(snapshot_idx, active_only=False)],
                self.receptor.n_atoms) > 0)).transpose()

            for pi, pj in contact_pair:
                p1 = self.pockets(snapshot_idx, active_only=False)[pi]
                p2 = self.pockets(snapshot_idx, active_only=False)[pj]
                if {p1.core_aux_minor, p2.core_aux_minor} in {{'core', 'aux'}, {'core'}, {'core', 'minor'},
                                                              {'aux', 'minor'}}:
                    pocket_vector1 = p1.lining_atoms_centroid - p1.centroid
                    pocket_vector2 = p2.lining_atoms_centroid - p2.centroid
                    if getCosAngleBetween(pocket_vector1, pocket_vector2) > 0:  # pocket vector facing inwards
                        pocket_graph.add_edge(p1, p2)

            core_aux_communities = [c for c in (networkx.connected_components(pocket_graph)) if
                                    len(c) > 1 and any([p.core_aux_minor == 'core' for p in c])]

            communities = []
            for community in core_aux_communities:
                linked_minor = set(chain.from_iterable([pocket_graph.neighbors(p) for p in community]))
                communities.append(linked_minor.union(community))
            self.communities[snapshot_idx] = communities

    def _gen_communities(self):

        # Todo convert to function and use multi core
        import networkx
        self.communities = {}
        self._pocket_network = {}
        for snapshot_idx in self.snapshots_indices:
            print(snapshot_idx)
            pocket_graph = networkx.Graph()
            pocket_graph.add_nodes_from(self.pockets(snapshot_idx, active_only=False))

            contact_pair = np.array(np.where(combination_intersection_count(
                [pocket.lining_atoms_idx for pocket in self.pockets(snapshot_idx, active_only=False)],
                self.receptor.n_atoms) > 0)).transpose()

            for pi, pj in contact_pair:
                p1 = self.pockets(snapshot_idx, active_only=False)[pi]
                p2 = self.pockets(snapshot_idx, active_only=False)[pj]
                if "core" in {p1.core_aux_minor, p2.core_aux_minor}:
                    pocket_vector1 = p1.lining_atoms_centroid - p1.centroid
                    pocket_vector2 = p2.lining_atoms_centroid - p2.centroid
                    if getCosAngleBetween(pocket_vector1, pocket_vector2) > 0:  # pocket vector facing inwards
                        pocket_graph.add_edge(p1, p2)

            core_aux_communities = [c for c in (networkx.connected_components(pocket_graph)) if
                                    len(c) > 1 and any([p.core_aux_minor == 'core' for p in c])]

            communities = []
            for community in core_aux_communities:
                linked_minor = set(chain.from_iterable([pocket_graph.neighbors(p) for p in community]))
                communities.append(linked_minor.union(community))
            self.communities[snapshot_idx] = communities
            self._pocket_network[snapshot_idx] = pocket_graph

    @property
    def snapshots_indices(self):
        return iter(range(self.n_frames))




    def screen_pockets(self):
        assert len(list(self.receptor._data.snapshots_idx)) == self.n_frames
        data = self.receptor._data

        if self.config.screen_by_face:
            interface_atom_idx = np.sort(self._get_face_atoms())
            for snapshot_idx in range(self.n_frames):
                for pocket in self.pockets(snapshot_idx):
                    atom_list = np.concatenate([pocket.lining_atoms_idx, interface_atom_idx])
                    if len(np.unique(atom_list)) == len(atom_list):
                        data[pocket.alpha_idx, 11] = 0

        if self.config.screen_by_space:
            assert self.config.min_space > 0
            for snapshot_idx in range(self.n_frames):
                for pocket in self.pockets(snapshot_idx):
                    if pocket.get_total_space < self.config.min_space:
                        pocket.deactivate()

        if self.config.min_num_alph > 0:
            for snapshot_idx in range(self.n_frames):
                for pocket in self.pockets(snapshot_idx):
                    if len(pocket.alpha_idx) <= self.config.min_num_alph:
                        pocket.deactivate()

    def _gen_d_pockets(self):
        """
        Generate d-pocket dictionary of list of indices
        :return: dict of d_pockets
        """
        if self._data is None:
            raise Exception('No Data available ')

        pockets_all = []
        for snapshot_idx in self.snapshots_indices:
            pockets_all.extend(self.pockets(snapshot_idx))

        lining_atom_indices = [pocket.lining_atoms_idx_redu for pocket in pockets_all]

        # calculate jaccard_diff_matrix

        intersection_matrix = combination_intersection_count(lining_atom_indices, self.receptor.n_atoms)

        union_matrix = combination_union_count(lining_atom_indices, self.receptor.n_atoms)

        jaccard_diff_matrix = 1 - intersection_matrix / union_matrix

        pocket_d_idx = list(fcluster(Z=linkage(squareform(jaccard_diff_matrix), method='average'),
                                     t=self.config.dpocket_cluster_cutoff,
                                     criterion='distance') - 1)

        d_pockets = {i: [] for i in range(max(pocket_d_idx) + 1)}

        for pocket_idx, dp_idx in zip(pockets_all, pocket_d_idx):
            d_pockets[dp_idx].append(pocket_idx)

        d_pockets = dict(list(enumerate(sorted(d_pockets.values(), reverse=True, key=len))))

        self._d_pockets = d_pockets

        return d_pockets

    def _gen_d_communities(self):
        if not self._d_pockets:
            self._gen_d_pockets()

        if not self._pocket_network:
            self._gen_communities()

        core_d_pockets = {}
        for i, pockets in self._d_pockets.items():
            for pocket in pockets:
                if pocket.core_aux_minor == "core":
                    core_d_pockets[i] = pockets
                    for p in pockets:
                        p._is_core_d = True
                    break

        # mark all pockets connected to members of d-pocket
        for snapshot_idx in self.snapshots_indices:
            network = self._pocket_network[snapshot_idx]
            for p in network.nodes():
                if p._is_core_d:
                    for cp in network[p]:
                        cp._connected = True

        for i, pockets in self._d_pockets.items():
            if i not in core_d_pockets:
                print(len([True for p in pockets if p._connected]))



    """
    Visualization methods
    """

    def view_snapshot(self, snapshot_idx: int = 0) -> object:
        self.view_receptor(snapshot_idx)
        self.view_binder(snapshot_idx)
        self.view_alphas(snapshot_idx, active_only=True)
        self.view_pocket_centers(snapshot_idx)
        self.view_pocket_surface(snapshot_idx)
        return self._view

    def view_binder(self, snapshot_idx=0):
        assert self._view
        self.binder_view = self._view.add_trajectory(self.binder.trajectory[snapshot_idx])

    def get_view(self):
        return self._view

    def view_receptor(self, snapshot_idx=0):
        try:
            import nglview as nv
        except:
            raise Exception('nglview is needed for jupyter notebook visualization')
        self._view = nv.show_mdtraj(self.receptor.trajectory[snapshot_idx], gui=True)
        self.receptor_view = self._view.component_0
        self.receptor_view.clear_representations()
        self.receptor_view.add_surface(selection='protein', opacity=0.8, color='white')

    # def show_pocket(self, snapshot_idx=0):
    #     self.pocket_view = self._view.add_trajectory(self.receptor.cluster(snapshot_idx).traj)
    #     self.pocket_view.clear_representations()
    #     self.pocket_view.add_representation(repr_type='ball+stick', selection='all', color='residueindex')

    def view_alphas(self, snapshot_idx=0, active_only=True):
        for pocket in self.pockets(snapshot_idx):
            color = self.config.color(idx=pocket._idx)
            if (not active_only) or (active_only and pocket.is_active):
                for alpha in pocket.alphas:
                    self._view.shape.add_buffer("sphere", position=list(alpha.xyz * 10), color=color, radius=[0.5])

    def view_pocket_centers(self, snapshot_idx=0, active_only=True):
        for pocket in self.pockets(snapshot_idx):
            color = self.config.color(idx=pocket._idx)
            if (not active_only) or (active_only and pocket.is_active):
                self._view.shape.add_buffer("sphere", position=list(pocket.centoid * 10), color=color, radius=[0.5])
                self._view.shape.add('text', list(pocket.centoid * 10), [0, 0, 0], 2.5, str(int(pocket._idx)))

    def view_pocket_surface(self, snapshot_idx=0):
        for pocket in self.pockets(snapshot_idx):
            self._view.add_surface(selection=list(pocket.lining_atoms_idx), opacity=1.0, color=pocket.color,
                                   surfaceType='sas')


if __name__ == '__main__':
    import mdtraj
    import sys

    test_receptor_path = sys.argv[1]
    test_binder_path = sys.argv[2]

    lig_traj = mdtraj.load(test_binder_path)
    prot_traj = mdtraj.load(test_receptor_path)

    complex = AS_Universe()
    complex.set_receptor(prot_traj)
    complex.set_binder(lig_traj)
