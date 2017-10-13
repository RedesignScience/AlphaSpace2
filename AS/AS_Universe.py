from AS.AS_Struct import *
from AS.AS_Config import AS_Config
from AS.AS_Cluster import *
import nglview as nv
import multiprocessing as mp
from mdtraj import shrake_rupley


# noinspection PyAttributeOutsideInit,PyAttributeOutsideInit,PyAttributeOutsideInit,PyAttributeOutsideInit,PyAttributeOutsideInit,PyAttributeOutsideInit
class AS_Session(object):
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
        if guess_receptor_binder and receptor is not None and binder is None:
            self.guess_receptor_binder(receptor, guess_by_order)
        else:
            self.set_receptor(receptor)
            self.set_binder(binder)
        self.others = None
        self.view = None
        self.n_frames = self.receptor.n_frames

    def __repr__(self):
        return "Receptor of {} residues {} atoms | Binder of {} residues {} atoms".format(self.receptor.n_residues,
                                                                                          self.receptor.n_atoms,
                                                                                          self.binder.n_residues,
                                                                                          self.binder.n_atoms)

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

    def pockets(self, snapshot_idx):
        return self.receptor.clusters[snapshot_idx].pockets

    def cluster(self, snapshot_idx=0):
        """
        return list of clusters
        :param snapshot_idx: int
        :return: object, AS_Cluster
        """
        return self.receptor.clusters[snapshot_idx]

    def clusters(self):
        """
        return list of clusters
        :return: list
        """
        return self.receptor.clusters

    def guess_receptor_binder(self, traj: object, by_order: bool = True) -> bool:
        """
        Divide receptor trajectory based on connectivity, set larger molecule as receptor.
        This process automatically removes water and other solvents
        :param traj: target trajectory
        :param by_order: bool, if False, guess by appearances in file
        :return: bool, if any macro molecule were found.
        """
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

    def set_binder(self, binder: object):
        """
        set binder (ligand) in session
        :param binder: object, trajectory
        :return:
        """
        self.binder = AS_Structure(binder, structure_type=1, parent=self)

    def set_receptor(self, receptor: object):
        """
        set receptor (protein) in session
        :param receptor: object, trajectory
        :return:
        """
        self.receptor = AS_Structure(receptor, structure_type=0, parent=self)

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
        self.receptor.generate_cluster(snapshot_idx=snapshot_idx)
        if self.binder:
            self.receptor.calculate_contact(binder=self.binder, snapshot_idx=snapshot_idx)

    def run_mp(self, cpu: int = 1):
        """
        run the AlphaSpace main program
        :param cpu: int, number of cpu you want to use, default use all
        """
        if cpu != 1:
            cpu = mp.cpu_count()
        pool = mp.Pool(cpu)
        pool.map(self.run, range(self.n_frames))

    def screen_by_ligand_contact(self, snapshot_idx: int = 0):
        self.receptor.clusters[snapshot_idx].screen_by_contact()

    def get_pockets(self, snapshot_idx: int = 0, force: bool = False) -> list:
        if len(self.receptor.clusters[snapshot_idx].pockets) == 0 or force:
            self.receptor.clusters[snapshot_idx]._build_pockets()
        return self.receptor.clusters[snapshot_idx].pockets

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
        self.view.component_2.add_representation(repr_type='label', lableType='residueindex')

    def show_binder(self, snapshot_idx=0):
        self.binder_view = self.view.add_trajectory(self.binder.trajectory[snapshot_idx])

    def show_pocket(self, snapshot_idx=0):
        self.pocket_view = self.view.add_trajectory(self.receptor.clusters[snapshot_idx].traj)
        self.pocket_view.clear_representations()
        self.pocket_view.add_representation(repr_type='ball+stick', selection='all', color='residueindex')

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

    def _get_face_pockets(self, snapshot_idx=0, interface_atom_index=None):
        """
        Get interface pockets
        :param snapshot_idx: int
        :return: list of AS_Pocket
        """
        pockets = self.get_pockets(snapshot_idx)
        if interface_atom_index is None:
            interface_atom_index = self._get_face_atoms(snapshot_idx)
        for pocket in pockets:
            pocket_interface_atom_index = interface_atom_index.intersection(pocket.get_lining_atoms())
            if len(pocket_interface_atom_index) / float(len(pocket.get_lining_atoms())) >= self.config.screen_face_perc:
                pocket.is_interface = True

        return [p for p in pockets if p.is_interface]






if __name__ == '__main__':
    import mdtraj

    test_binder_path = '/Users/haotian/Dropbox/pycharm_project/AlphaSpace/Test_system/lig.pdb'
    test_receptor_path = '/Users/haotian/Dropbox/pycharm_project/AlphaSpace/Test_system/prot.pdb'

    lig_traj = mdtraj.load(test_binder_path)
    prot_traj = mdtraj.load(test_receptor_path)

    complex = AS_Session()
    complex.set_receptor(prot_traj)
    complex.set_binder(lig_traj)
