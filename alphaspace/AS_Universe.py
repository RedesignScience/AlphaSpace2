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

from itertools import chain

import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from mdtraj.geometry.sasa import _ATOMIC_RADII
import mdtraj

from .AS_Cluster import AS_D_Pocket
from .AS_Config import AS_Config
from .AS_Funct import _tessellation_mp, getCosAngleBetween, combination_intersection_count_mp, \
    combination_intersection_count, combination_union_count
from .AS_Struct import AS_Structure


# noinspection PyAttributeOutsideInit,PyAttributeOutsideInit,PyAttributeOutsideInit,PyAttributeOutsideInit,PyAttributeOutsideInit,PyAttributeOutsideInit
class AS_Universe(object):
    def __init__(self, receptor=None, binder=None, guess_receptor_binder=True, guess_by_order=True, config=None,
                 tag="", keepH=False):
        """
        Container for an AlphaSpace session, have child container receptor and binder
        :param receptor: object
        :param binder: object
        :param guess_by_order: book
        :param guess_receptor_binder: bool, guess is based on molecule size
        :param config: object,AS_config
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
        self._view = None

        self._d_pockets = {}
        self._pocket_network = {}

        self.tag = tag

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
        TODO
        calculate the d pockets and give a iterator os AS_D_Pocket
        :return: object, AS_D_Pocket
        """
        if not self._d_pockets:
            self._gen_d_pockets()
        for p in self._d_pockets:
            yield self._d_pockets[p]

    def d_pocket(self, i) -> AS_D_Pocket:
        """
        todo
        :param i:
        :return:
        """
        if not self._d_pockets:
            self._d_pockets = self.receptor._gen_d_pockets()
        return self._d_pockets[i]

    # def _is_processed(self, snapshot_idx: int) -> bool:
    #     if snapshot_idx in self.receptor._clusters:
    #         return True
    #     else:
    #         return False

    def pockets(self, snapshot_idx: int = 0, active_only: bool = False) -> list:
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
            x = self.binder.trajectory + structure
            self.binder.trajectory = x
        else:
            self.binder = AS_Structure(structure, structure_type=1, parent=self)

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
            x = self.receptor.trajectory + structure
            self.receptor.trajectory = x
        else:
            self.receptor = AS_Structure(structure, structure_type=0, parent=self)

    def run_alphaspace(self):
        self.run_alphaspace_mp(cpu=1)

    def run_alphaspace_mp(self, cpu=None):
        _tessellation_mp(self, cpu=cpu)

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
        self._communities = {}
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
            self._communities[snapshot_idx] = communities

    def _gen_communities(self):

        # Todo convert to function and use multi core
        import networkx
        self._communities = {}
        self._pocket_network = {}
        for snapshot_idx in self.snapshots_indices:
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
            self._communities[snapshot_idx] = communities
            self._pocket_network[snapshot_idx] = pocket_graph

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
        # check if alphaspace is run

        if self.data is None:
            raise Exception('No Data available ')

        # list all pockets
        pockets_all = []
        for snapshot_idx in self.snapshots_indices:
            pockets_all.extend(self.pockets(snapshot_idx))

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

    def calculate_vina_score(self, snapshot_idx=0):
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
        for pocket in self.pockets(active_only=True):
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
            beta.prb_element = prb_element
            beta.vina_score = prb_score[i]

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
