import nglview as nv
import numpy as np
# import multiprocessing as mp
import pathos.multiprocessing as mp
from mdtraj import shrake_rupley
from AS_Config import AS_Config
from AS_Struct import AS_Structure
from AS_Cluster import AS_Cluster,AS_D_Pocket
from AS_Funct import _tessellation

from asyncio import get_event_loop, wait, ensure_future


# noinspection PyAttributeOutsideInit,PyAttributeOutsideInit,PyAttributeOutsideInit,PyAttributeOutsideInit,PyAttributeOutsideInit,PyAttributeOutsideInit
class AS_Universe(object):
    def __init__(self, receptor=None, binder=None, guess_receptor_binder=True, guess_by_order=True, config=None):
        """
        Container for an AlphaSpace session, have child container receptor and binder
        :param receptor: object
        :param binder: object
        :param guess_by_order: book
        :param guess_receptor_binder: bool, guess is based on molecule size
        :param config: object,AS_config
        """

        if config is not None:
            self.config = config
        else:
            self.config = AS_Config()

        self.set_receptor(receptor)
        self.set_binder(binder)

        if guess_receptor_binder and receptor:
            self.guess_receptor_binder(receptor, guess_by_order)

        self.others = None
        self.view = None

        self._d_pockets = {}



    def __repr__(self):
        return "Receptor of {} residues {} atoms | Binder of {} residues {} atoms".format(self.receptor.n_residues,
                                                                                          self.receptor.n_atoms,
                                                                                          self.binder.n_residues,
                                                                                          self.binder.n_atoms)
    @property
    def frames(self):
        for i in range(self.n_frames):
            yield i


    @property
    def clusters(self):
        """
        return list of clusters
        :return: list
        """
        return self.receptor.clusters

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

    def d_pocket(self,i) -> AS_D_Pocket:
        if not self._d_pockets:
            self._d_pockets = self.receptor._gen_d_pockets()
        return self._d_pockets[i]

    def _is_processed(self, snapshot_idx: int) -> bool:
        if snapshot_idx in self.receptor._clusters:
            return True
        else:
            return False

    def pockets(self, snapshot_idx):
        return self.receptor._clusters[snapshot_idx].pockets

    def cluster(self, snapshot_idx: int = 0) -> object:
        """
        return list of clusters
        :param snapshot_idx: int
        :return: object, AS_Cluster
        """
        return self.receptor._clusters[snapshot_idx]

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
            self.set_receptor(traj.atom_slice([atom.index for atom in molecule_list[0]]))
            self.set_binder(traj.atom_slice([atom.index for atom in molecule_list[1]]))
            return True
        elif len(molecule_list) == 1:
            self.set_receptor(traj.atom_slice([atom.index for atom in molecule_list[0]]))
            self.binder = None
            return True
        else:
            return False

    def set_binder(self, structure: object, append=False):
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

    def set_receptor(self, structure: object, append=False):
        """
        set receptor (protein) in session
        :param structure: object, trajectory
        :param append: Bool, if the new binder should be appended to the preview one, default overwritten.
        :return:
        """
        if structure is None:
            self.receptor = None
            return

        if append and (self.receptor is not None):
            x = self.receptor.trajectory + structure
            self.receptor.trajectory = x

        else:
            self.receptor = AS_Structure(structure, structure_type=0, parent=self)

    def set_others(self, others: object):
        """
        set other molecules: water etc
        :param others: object, trajectory
        :return:
        """
        self.others = AS_Structure(others, structure_type=0, parent=self)

    def run(self, snapshot_idx=0):
        """
        Private method, please use run
        """
        cluster = self.receptor.generate_cluster(snapshot_idx=snapshot_idx)
        cluster._tessellation(self.config)
        if self.binder is not None:
            self.receptor.calculate_contact(binder=self.binder, snapshot_idx=snapshot_idx)
        return self.cluster(snapshot_idx)

    def run_all(self):
        """
        run the AlphaSpace main program
        :param cpu: int, number of cpu you want to use, default use all
        """

        pool = mp.Pool()
        results = pool.map(self.run, range(self.n_frames))








    def screen_by_ligand_contact(self, snapshot_idx: int = 0):
        self.receptor.clusters[snapshot_idx].screen_by_contact()

    def get_pockets(self, snapshot_idx: int = 0) -> list:

        return self.pockets(snapshot_idx)

    def _get_face_atoms(self, snapshot_idx=0):
        """
        Calculate the snapshot interface atom.
        The interface atom is defined as whose ASA is reduced with introduction of ligand.
        :param snapshot_idx: int
        :return: numpy.array
        """
        receptor_snapshot = self.receptor.traj[snapshot_idx]

        complex_snapshot = receptor_snapshot.stack(self.binder.traj[snapshot_idx])

        receptor_snapshot_sasa = shrake_rupley(receptor_snapshot)[0]

        complex_snapshot_sasa = shrake_rupley(complex_snapshot)[0]

        sasa_diff = receptor_snapshot_sasa - complex_snapshot_sasa[:len(receptor_snapshot_sasa)]

        interface_atom_index = np.where(sasa_diff > 0)[0]

        return set(interface_atom_index)

    def screen_pockets(self):
        if self.receptor.n_clusters < 1:
            raise Exception("Must process all frames before screening pockets.")
        if self.config.screen_by_lig_cntct:
            for snapshot_idx in range(self.n_frames):
                contact = self.cluster(snapshot_idx)._get_contact_list()
                self.cluster(snapshot_idx)._slice(np.where(contact != 0)[0])

        if self.config.screen_by_face and not self.config.screen_by_lig_cntct:
            for snapshot_idx in range(self.n_frames):
                face_pocket_alpha_index = set()
                face_atoms_index = self._get_face_atoms(snapshot_idx)
                for pocket in self.cluster(snapshot_idx).pockets:
                    lining_atoms = self.cluster(snapshot_idx)._get_lining_atoms(pocket.get_alpha_index())
                    if len(face_atoms_index.intersection(lining_atoms)) > 0:
                        face_pocket_alpha_index = face_pocket_alpha_index.union(set(pocket.get_alpha_index()))
                self.cluster(snapshot_idx)._slice(face_pocket_alpha_index)

        if self.config.screen_by_score:
            for cluster in self.receptor.clusters:
                face_pocket_alpha_index = set()
                for pocket in cluster.pockets:
                    if pocket.get_total_score > self.config.min_score:
                        face_pocket_alpha_index.union(set(pocket.get_alpha_index))
                cluster._slice(face_pocket_alpha_index)

        if self.config.screen_by_res:
            for cluster in self.receptor.clusters:
                face_pocket_alpha_index = set()
                # TODO Finish this section
                for pocket in cluster.pockets:
                    if pocket.get_lining_residues:
                        face_pocket_alpha_index.union(set(pocket.get_alpha_index))
                cluster._slice(face_pocket_alpha_index)

        self._d_pockets = None

    """
    Visualization methods
    """

    def view_snapshot(self, snapshot_idx: int = 0) -> object:
        self.show_receptor(snapshot_idx)
        self.show_binder(snapshot_idx)
        self.show_pocket(snapshot_idx)
        return self.view

    def show_receptor(self, snapshot_idx=0):
        self.view = nv.show_mdtraj(self.receptor.trajectory[snapshot_idx], gui=True)
        self.receptor_view = self.view.component_0
        self.receptor_view.clear_representations()
        self.receptor_view.add_surface(selection='protein', opacity=1, color='white')

    def show_pocket_label(self):
        self.view.component_2.add_representation(repr_type='label', lableType='residueindex', color='residueindex')

    def show_binder(self, snapshot_idx=0):
        self.binder_view = self.view.add_trajectory(self.binder.trajectory[snapshot_idx])

    def show_pocket(self, snapshot_idx=0):
        self.pocket_view = self.view.add_trajectory(self.receptor.cluster(snapshot_idx).traj)
        self.pocket_view.clear_representations()
        self.pocket_view.add_representation(repr_type='ball+stick', selection='all', color='residueindex')




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
