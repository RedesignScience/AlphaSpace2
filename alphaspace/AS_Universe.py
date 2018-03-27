"""
The main container AS_Universe
==============================

When creating a new alphaspace instance, all the data are stored in an object being referred to as universe.
And nearly all functions can be done by calling methods in the AS_Universe.

Initialization
--------------

You can initialize the universe object by calling it and passing receptor and binder objects.

>>> from alphaspace import *
>>> import mdtraj
>>> import sys
>>>
>>> universe = AS_Universe(receptor=mdtraj.load(sys.argv[1),binder=mdtraj.load(sys.argv[2))


AS_Universe can intelligently construct molecules from your input structure and guess which one is the receptor or
binder. This is implemented that by default, the largest molecule (one with most atoms) is considered to be receptor,
the second largest is considered binder, if that molecule is not a solvent or ion.

If you wish to guess the receptor and binder based on their appearance in the file:

>>> universe = AS_Universe(receptor=mdtraj.load(sys.argv[1),binder=mdtraj.load(sys.argv[2),
>>>                                         guess_by_order=True )

If you wish to define the receptor and binder yourself, or the guessing algorithm does not work properly,
you can load them separately.
You can set the receptor and binder after you have initialized the universe.

>>> universe = AS_Universe()
>>> universe.set_receptor(mdtraj.load(sys.argv[1))
>>> universe.set_binder(mdtraj.load(sys.argv[2])

You can checkout more on how to load and select molecular object in the receptor_binder_selection.py


"""

import numpy as np

from .AS_Cluster import AS_D_Pocket
from .AS_Config import AS_Config
from .AS_Funct import _tessellation_mp, getCosAngleBetween, combination_intersection_count
from .AS_Struct import AS_Structure

import nglview as nv

import networkx as nx
from itertools import combinations, combinations_with_replacement
from .AS_Funct import is_pocket_connected


# noinspection PyAttributeOutsideInit,PyAttributeOutsideInit,PyAttributeOutsideInit,PyAttributeOutsideInit,PyAttributeOutsideInit,PyAttributeOutsideInit
class AS_Universe(object):
    def __init__(self, receptor=None, binder=None, guess_receptor_binder=True, guess_by_order=True, config=None,
                 label="", keepH=False):
        """

        Parameters
        ----------
        receptor
        binder
        guess_receptor_binder
        guess_by_order
        config
        label
        keepH

        Returns
        -------
        universe : AS_Universe
        """

        self.config = config if config is not None else AS_Config()

        self.receptor = None
        self.binder = None
        self.set_receptor(receptor, keepH)
        self.set_binder(binder, keepH)

        if guess_receptor_binder and receptor and not binder:
            if self.guess_receptor_binder(receptor, guess_by_order):
                pass
            else:
                raise Exception('No binder detected')

        self._pocket_list = {}

        self.others = None
        self._ngl_view = None

        self._d_pockets = {}
        self._pocket_network = {}

        self.label = label

        self._ngl_added_component = []

    def __repr__(self):
        rec_res = self.receptor.n_residues if self.receptor else 0
        rec_atm = self.receptor.n_atoms if self.receptor else 0
        bind_res = self.binder.n_residues if self.binder else 0
        bind_atm = self.binder.n_atoms if self.binder else 0

        return "Receptor of {} residues {} atoms | Binder of {} residues {} atoms".format(rec_res,
                                                                                          rec_atm,
                                                                                          bind_res,
                                                                      bind_atm)





    @property
    def data(self):
        return self.receptor._data

    @property
    def frames(self):
        for i in range(self.n_frames):
            yield i

    @property
    def frame_indices(self):
        return self.frames

    @property
    def n_frames(self):
        if self.receptor is None:
            return 0
        else:
            return self.receptor.trajectory.n_frames

    @property
    def n_atoms(self):
        """
        return the total number of atoms in receptor and binders

        Returns
        -------

        int

        """
        """
        
        :return: 
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
            self._gen_d_pockets()
        for p in self._d_pockets:
            yield self._d_pockets[p]

    def d_pocket(self, i) -> AS_D_Pocket:
        """

        Parameters
        ----------
        i

        Returns
        -------

        """
        if not self._d_pockets:
            self._d_pockets = self.receptor._gen_d_pockets()
        return self._d_pockets[i]

    def pockets(self, snapshot_idx: int = 0, active_only: bool = False) -> list:
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
            else:
                continue

        if not by_order:
            molecule_list.sort(key=len, reverse=True)

        if len(molecule_list) > 1:
            self.set_receptor(traj.atom_slice(np.array([atom.index for atom in molecule_list[0]], dtype=int)))
            self.set_binder(traj.atom_slice(np.array([atom.index for atom in molecule_list[1]], dtype=int)))
            return True
        elif len(molecule_list) == 1:
            self.set_receptor(traj.atom_slice([atom.index for atom in molecule_list[0]]))
            self.binder = None
            return False
        else:
            return False

    def set_binder(self, structure, append=False, keepH=False):

        """
        set binder (ligand) in session
        :param structure: object, trajectory
        :param append: Bool, if the new binder should be appended to the preview one, default overwritten.
        :return:
        """

        from mdtraj.core import element
        if structure is None:
            self.binder = None
            return

        if not keepH:
            non_h_idx = [a.index for a in structure.topology.atoms if a.element != element.hydrogen]
            structure.atom_slice(non_h_idx, inplace=True)

        if append and (self.binder is not None):
            x = self.binder.trajectory.join(structure, check_topology=True)
            self.binder.trajectory = x
        else:
            self.binder = AS_Structure(structure, structure_type=1, parent=self)

    def remove_h(self):
        from mdtraj.core import element

        if self.receptor:
            non_h_idx = [a.index for a in self.receptor.topology.atoms if a.element != element.hydrogen]
            self.receptor.traj.atom_slice(non_h_idx, inplace=True)

        if self.binder:
            non_h_idx = [a.index for a in self.binder.topology.atoms if a.element != element.hydrogen]
            self.binder.traj.atom_slice(non_h_idx, inplace=True)

    def set_receptor(self, structure, append=False, keepH=False):
        """
        set receptor (protein) in session
        :param structure: trajectory
        :param append: Bool, if the new binder should be appended to the preview one, default overwritten.
        :return:
        """

        if structure is None:
            self.receptor = None
            return False

        if not keepH:
            non_h_idx = structure.topology.select_atom_indices(selection='heavy')
            structure.atom_slice(non_h_idx, inplace=True)

        if append and (self.receptor is not None):
            x = self.receptor.trajectory.join(structure, check_topology=True)
            self.receptor.trajectory = x
        else:
            self.receptor = AS_Structure(structure, structure_type=0, parent=self)

    def run_alphaspace(self, frame_range=None):
        self.run_alphaspace_mp(cpu=1, frame_range=frame_range)

    def run_alphaspace_mp(self, frame_range=None, cpu=None, ):
        _tessellation_mp(self, cpu=cpu, frame_range=frame_range)

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



    def _gen_communities_legacy(self):
        import networkx
        self._communities = {}
        for snapshot_idx in range(self.n_frames):
            pocket_graph = networkx.Graph()
            pocket_graph.add_nodes_from(self.pockets(snapshot_idx, active_only=False))

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

            communities = []
            for cp in pocket_graph.nodes():
                community = [cp].extend(pocket_graph[cp])
                communities.append(community)

            # for community in core_aux_communities:
            #     linked_minor = set(chain.from_iterable([pocket_graph.neighbors(p) for p in community]))
            #     communities.append(linked_minor.union(community))
            self._communities[snapshot_idx] = communities

    def _connect_pockets(self, snapshot_idx=0):
        """
        Generate a networkx graph with nodes from pockets in the given snapshot.
        Parameters
        ----------
        snapshot_idx
        int

        Returns
        -------
        pocket_graph
            undirected graph

        """
        pocket_graph = nx.Graph()
        pocket_graph.add_nodes_from(self.pockets(snapshot_idx, False))
        for p1, p2 in combinations(pocket_graph.nodes(), 2):
            if is_pocket_connected(p1, p2):
                pocket_graph.add_edge(p1, p2)

        return pocket_graph

    def _gen_community(self, core_cutoff=100):
        """

        Parameters
        ----------
        core_cutoff
        int : space for core cutoff

        """
        from .AS_Cluster import AS_community

        self._communities = {}
        for snapshot_idx in self.snapshots_indices:
            pockets = sorted(self.pockets(snapshot_idx), reverse=True, key=lambda p: p.space)

            communities = []

            if pockets[0].space < 100:
                cores = [pockets[0]]
            else:
                cores = [pocket for pocket in pockets if pocket.space > core_cutoff]

            if len(cores) > 1:
                core_graph = nx.Graph()
                core_graph.add_nodes_from(cores)
                for p1, p2 in combinations(core_graph.nodes(), 2):
                    if is_pocket_connected(p1, p2):
                        core_graph.add_edge(p1, p2)

                for com in nx.connected_components(core_graph):
                    communities.append(AS_community(list(com)))

            for p in pockets[1:]:
                if p.space < core_cutoff:
                    for c in communities:
                        if c.connected_core(p):
                            c.aux.append(p)

            self._communities[snapshot_idx] = communities

    def communities(self, snapshot_idx=0):
        """
        Get communities from a snapshot
        Parameters
        ----------
        snapshot_idx

        Returns
        -------
        communities : list
        """
        return self._communities[snapshot_idx]

    @property
    def snapshots_indices(self):
        return iter(range(self.n_frames))

    def screen_pockets(self):
        assert len(list(self.data.snapshots_idx)) == self.n_frames

        if self.config.screen_by_face:
            interface_atom_idx = np.sort(self._get_face_atoms())
            for snapshot_idx in range(self.n_frames):
                for pocket in self.pockets(snapshot_idx):
                    atom_list = np.concatenate([pocket.lining_atoms_idx, interface_atom_idx])
                    if len(np.unique(atom_list)) == len(atom_list):
                        self.data[pocket.alpha_idx, 11] = 0

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
        Generate d pockets

        Returns
        -------

        d_pocket: dict
            a dict of d pockets, each is a list of child pocket indices

        """
        from alphaspace.AS_Funct import cluster_by_overlap

        if self.data is None:
            raise Exception('No Data available ')

        # list all pockets
        pockets_all = []
        lining_atom_indices = []
        for snapshot_idx in self.snapshots_indices:
            for pocket in self.pockets(snapshot_idx):
                if pocket.space > 20:
                    pockets_all.append(pocket)
                    lining_atom_indices.append(pocket.lining_atoms_idx)
        print("Clustering pockets")

        # extract lining atom into list
        lining_atom_indices = [pocket.lining_atoms_idx for pocket in pockets_all]

        d_pocket_p_idx = cluster_by_overlap(lining_atom_indices, self.receptor.n_atoms,
                                            self.config.dpocket_cluster_cutoff)

        d_pockets = dict(list(enumerate(sorted(d_pocket_p_idx.values(), reverse=True, key=len))))
        # fill pocket list
        for key, idx in d_pockets.items():
            d_pockets[key] = [pockets_all[i] for i in idx]

        self._d_pockets = d_pockets
        return d_pockets

    def _gen_d_pockets_iter(self, sample_frames=20, sample_ratio=1, pocket_space_cutoff=20):
        """

        Leader follower iterative generation of dpockets to save memory and improve speed.


        Parameters
        ----------
        sample_frames : int
            how many frames to sample for leaders
        sample_ratio : float
            the sampling fraction of pockets in leader d_pockets.
            1 for keeping all, 10 for picking 1/10 of pockets in each dpocket
        pocket_space_cutoff : int
            cutoff for pocket because small pockets are not considered.

        Returns
        -------
        d_pocket : dict
            a dictionary of lists
            each list contains pockets in the d_pocket.

        """
        from .AS_Funct import cluster_by_overlap, _prune_dpockets
        from collections import defaultdict

        # sample 20 snapshots from the universe
        sampled_snapshot_idx = np.random.choice(self.n_frames, sample_frames, False)
        leader_pockets = []
        for ss_idx in sampled_snapshot_idx:
            for pocket in self.pockets(ss_idx, active_only=False):
                if pocket.space > pocket_space_cutoff:
                    leader_pockets.append(pocket)

        # first cluster these frames' pockets into d_pockets

        d_pocket_p_idx = cluster_by_overlap([pocket.lining_atoms_idx for pocket in leader_pockets],
                                            self.receptor.n_atoms,
                                            self.config.dpocket_cluster_cutoff)

        # Create leader pockets dictionary
        leader_dpocket_dict = {}
        for dp_idx, p_idx in d_pocket_p_idx.items():
            leader_dpocket_dict[dp_idx] = [leader_pockets[i] for i in p_idx]

        # sample and generate a pruned pocket list with labels of leader index.

        leader_pockets, leader_label = _prune_dpockets(leader_dpocket_dict, sample_ratio)

        # create array for all pockets' lining atoms
        leader_pockets_atom_array = np.array([pocket.lining_atoms_vector for pocket in leader_pockets])

        follower_dpocket_dict = defaultdict(lambda: [])

        # iterate through snapshots

        for ss_idx in self.snapshots_indices:
            new_pocket_d_idx = []
            current_pockets = (p for p in self.pockets(snapshot_idx=ss_idx, active_only=False) if
                               p.space > pocket_space_cutoff)

            for pocket in current_pockets:
                pocket_lining_atom_array = pocket.lining_atoms_vector
                overlap = np.dot(leader_pockets_atom_array, pocket_lining_atom_array)
                union = np.count_nonzero(leader_pockets_atom_array + pocket_lining_atom_array, axis=1)
                tanimoto = overlap / union

                max_idx = np.argmax(tanimoto)

                new_pocket_d_idx.append(leader_label[int(max_idx)])

                follower_dpocket_dict[leader_label[int(max_idx)]].append(pocket)
            print("clustering dpocket in {}".format(ss_idx))

        self._d_pockets = follower_dpocket_dict

        return follower_dpocket_dict

    def _connect_d_pockets(self):
        """
        This graph is build on D-pockets detected in the snapshots.

        Exchange_overlap = Union_overlap - (Time Overlap)

        Union_overlap is overlap count for all lining atom set in dpocket

        Time Overlap is the overlap of each pockets' lining atom in each snapshot.

        Returns
        -------
        dpocket_graph
            networkx graph for d-pockets

        """
        # calculate the lining atom matrix

        d_pocket_lining_atom_matrices = []

        for dp in self.d_pockets:
            matrix = np.zeros((self.n_frames, self.receptor.n_atoms))
            for p in dp:
                np.put(matrix[p.snapshot_idx], p.lining_atoms_idx, 1)
            d_pocket_lining_atom_matrices.append(matrix)

        num_dp = len(d_pocket_lining_atom_matrices)

        # calculate the union overlap
        d_pocket_lining_atom_unions = [
            np.sum(m, axis=0, dtype=bool) for m in d_pocket_lining_atom_matrices
        ]

        union_overlap = np.zeros((num_dp, num_dp), dtype=float)
        time_overlap = np.zeros((num_dp, num_dp), dtype=float)
        for i, j in combinations(range(num_dp), 2):
            # calculate the union overlap
            union_overlap_vector = np.logical_and(d_pocket_lining_atom_unions[i], d_pocket_lining_atom_unions[j])
            union_overlap[i, j] = union_overlap[j, i] = np.count_nonzero(union_overlap_vector)
            # calculate the time overlap
            time_overlap_vector = np.logical_and(d_pocket_lining_atom_matrices[i], d_pocket_lining_atom_matrices[j])
            time_overlap[i, j] = time_overlap[j, i] = np.average(np.sum(time_overlap_vector, axis=1, dtype=float))

        return union_overlap, time_overlap

    def _gen_d_communities(self):
        if not self._d_pockets:
            self._gen_d_pockets()

        if not self._pocket_network:
            self._gen_community()

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

        # for i, pockets in self._d_pockets.items():
        #     if i not in core_d_pockets:
        #         print(len([True for p in pockets if p._connected]))

    def load(self, as_file):
        import _pickle as pickle

        with open(as_file, 'rb') as handle:
            u = pickle.load(handle)

        self.__dict__.update(u.__dict__)

    def dump(self, as_file):

        import _pickle as pickle

        with open(as_file, 'wb') as handle:
            pickle.dump(self, handle)

    def set_pdbqt(self, pdbqt_file):
        """
        Load pdbqt file of the receptor protein.

        Parameters
        ----------
        pdbqt_file : str
            path of pdbqt file
        """
        from .AS_Vina import pre_process_pdbqt
        self.pdbqt_prot_coord, self.prot_types, self.hp_type, self.acc_type, self.don_type = pre_process_pdbqt(
            pdbqt_file, truncation_length=self.receptor.n_atoms)

    def calculate_vina_score(self, snapshot_idx=0, active_only=False):
        """
        Calculate the vina score for all beta atoms in the given snapshot.
        This action requires:

        1. Main tessellation was performed on the given snapshot, so beta atoms can be generated.
        2. Protein atom types must be specified by assigning pdbqt file in `.load_pdbqt()`

        After finishing this calculation, you can access the score through beta atom or pocket via method:
        'beta.score' or 'pocket.score'

        See Also
        --------
        load_pdbqt : load atom type from pdbqt file.


        Parameters
        ----------
        snapshot_idx : int
            index of snapshot you wish to calculate
        """
        from .AS_Vina import get_probe_score

        betas = []
        probe_coords = []
        for pocket in self.pockets(snapshot_idx, active_only):
            for beta in pocket.betas:
                betas.append(beta)
                probe_coords.append(beta.centroid)

        prb_dict = get_probe_score(self.receptor.traj.xyz[snapshot_idx] * 10, self.prot_types, self.hp_type,
                                   self.acc_type, self.don_type,
                                   probe_coords=np.array(probe_coords) * 10)

        prb_score = []
        prb_element = []

        for i in prb_dict:
            prb_element.append(i)
            prb_score.append(prb_dict[i])

        prb_score = np.array(prb_score).transpose((1, 0, 2))

        for i, beta in enumerate(betas):
            # beta.prb_element = prb_element
            beta._vina_score = prb_score[i]

    """
    Visualization methods
    """

    def draw(self, snapshot_idx: int = 0, show_binder=True, surface_opacity=1.0):
        """
        Draw and view the current snapshot in jupyter notebook.

        Parameters
        ----------
        surface_opacity
        snapshot_idx : int
        show_binder : bool

        Returns
        -------
        view
            nglview widget object
        """

        self._ngl_view = nv.show_mdtraj(self.receptor.traj)

        if surface_opacity > 0:
            self._ngl_view.add_surface(selection='protein', opacity=surface_opacity, color='white')

        if self.binder and show_binder:
            self._ngl_view.add_trajectory(self.binder.traj)

        self._ngl_view.frame = snapshot_idx

        self._ngl_added_component = []
        return self._ngl_view

    @property
    def view(self):
        """
        Get the nglview object

        Returns
        -------

        """
        if self._ngl_view:

            return self._ngl_view

        else:
            raise Exception("need to generate view first with .draw method")

    def _draw_as_sphere(self, item, color, radius=None, opacity=1.0):
        """

        Parameters
        ----------
        item alpha_atom or xyz coordinates in Angstrom
        color
        radius
        opacity

        Returns
        -------

        """
        if not self._ngl_added_component:
            self._ngl_added_component = list(range(self._ngl_view.n_components))

        if type(item) == list:
            self._ngl_view.shape.add_buffer("sphere", position=item, color=color, radius=[radius])
        else:
            _radius = radius if radius is not None else item._ngl_radius
            self._ngl_view.shape.add_buffer("sphere", position=list(item.centroid * 10), color=color,
                                            radius=[_radius])
        _ngl_component_idx = self._ngl_added_component[-1] + 1
        self._ngl_added_component.append(_ngl_component_idx)
        self._ngl_view._remote_call('updateRepresentationForComponent', target='Widget',
                                    args=[0, _ngl_component_idx], kwargs={'opacity': opacity})
        return _ngl_component_idx

    def _draw_cylinder(self, position1, position2, color, radius, opacity=1.0):

        if not self._ngl_added_component:
            self._ngl_added_component = list(range(self._ngl_view.n_components))

        _ngl_component_idx = self._ngl_added_component[-1] + 1
        self._ngl_added_component.append(_ngl_component_idx)

        self._ngl_view.shape.add_buffer("cylinder", position1=position1, position2=position2, color=color,
                                        radius=radius)
        self._ngl_view._remote_call('updateRepresentationForComponent', target='Widget', args=[0, _ngl_component_idx],
                                    kwargs={'opacity': opacity})
        return _ngl_component_idx


    def draw_pocket_graph(self, snapshot_idx=0, active_only=True):
        pocket_graph = self._connect_pockets(snapshot_idx)
        for p1, p2 in pocket_graph.edges:
            if active_only and (p1.is_active and p2.is_active):
                self._draw_cylinder(position1=list(p1.centroid * 10), position2=list(p2.centroid * 10),
                                    color=[0, 0, 0], radius=[0.1])
            else:
                self._draw_cylinder(position1=list(p1.centroid * 10), position2=list(p2.centroid * 10),
                                    color=[0, 0, 0], radius=[0.1])

    def _hide_item(self, item):
        if item is int:
            self._ngl_view._remote_call('removeRepresentation', target='Widget', args=[item._ngl_component_idx, 0])
        else:
            try:
                self._ngl_view._remote_call('removeRepresentation', target='Widget', args=[item, 0])
            except:
                raise TypeError('Item is of the wrong type')

    def view_alphas(self, snapshot_idx=0, active_only=True, opacity=1, show_noise=False):
        for pocket in self.pockets(snapshot_idx, active_only):
            if not show_noise and pocket.is_noise:
                continue
            color = pocket.color
            for alpha in pocket.alphas:
                self._draw_as_sphere(alpha, color=color, opacity=opacity)

    def view_betas(self, snapshot_idx=0, active_only=True, opacity=1, show_noise=False):

        for pocket in self.pockets(snapshot_idx, active_only):
            if not show_noise and pocket.is_noise:
                continue
            color = pocket.color
            for beta in pocket.betas:
                self._draw_as_sphere(beta, color=color, opacity=opacity)

    def view_pocket_centers(self, snapshot_idx=0, active_only=True, opacity=1, show_noise=False):

        for pocket in self.pockets(snapshot_idx, active_only):
            if not show_noise and pocket.is_noise:
                continue
            color = pocket.color
            self._draw_as_sphere(pocket, color=color, opacity=opacity)

    def clear_view(self):
        self._ngl_view.clear_representations()


def load(traj, top=None, label=None):
    """
    Load a file into alphaspace and return the AS_Universe object

    Parameters
    ----------

    file : str
        file path

    Returns
    -------
    universe : AS_Universe

    """
    import mdtraj as md

    traj = md.load(traj, top=top)

    return AS_Universe(receptor=traj, guess_receptor_binder=True, guess_by_order=False, label=label)


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
