import numpy as np
from scipy.spatial import Voronoi,Delaunay
from scipy.cluster.hierarchy import linkage,fcluster
from AS.AS_Funct import *
from mdtraj.core.topology import Topology
from mdtraj.core.trajectory import Trajectory


# noinspection PyAttributeOutsideInit,PyAttributeOutsideInit,PyAttributeOutsideInit,PyAttributeOutsideInit,PyAttributeOutsideInit,PyAttributeOutsideInit,PyTypeChecker
class AS_Cluster(object):
    def __init__(self,receptor,config,snapshot_idx = 0,parent = None):
        """
        Container for alpha, beta, gamma atoms
        :param receptor_snapshot: coordinates of parent receptor
        :param config: configuration file
        """
        self.config = config
        self.receptor_traj = receptor
        self.receptor_top = receptor.top
        self.snapshot_idx = snapshot_idx
        self.parent = None
        self.top = Topology()
        # create empty initial chain and residue as atom containers
        self.top.add_chain()
        self.top.add_residue(name='ASC',chain=self.top.chain(0))
        xyz = np.zeros([1,0,3])
        self.traj = Trajectory(xyz,self.top)
        self.atom_polarity = np.array(
            [(str(atom.element) in ['nitrogen', 'oxygen', 'sulfur']) for atom in self.receptor_top.atoms])
        self.tessellation(self.config)
        self.pockets = []

    def __repr__(self):
        return "Alpha Atom cluster of #{} frame".format(self.snapshot_idx)

    def tessellation(self, config):
        """
        perform tessellation in order to generate the cluster of alpha atoms.
        :param config: object
        """
        # Generate Raw Tessellation vertices
        raw_alpha_lining = Delaunay(self.receptor_traj.xyz[0]).simplices
        raw_alpha_xyz = Voronoi(self.receptor_traj.xyz[0]).vertices
        raw_lining_xyz = np.take(self.receptor_traj.xyz[0],raw_alpha_lining[:,0].flatten(),axis=0)

        # Calculate Raw alpha sphere radii
        raw_alpha_sphere_radii = np.linalg.norm(raw_lining_xyz - raw_alpha_xyz,axis=1)

        # Filter the data based on radii
        filtered_idx = np.where(np.logical_and(config.min_r / 10.0 <= raw_alpha_sphere_radii,raw_alpha_sphere_radii <= config.max_r / 10.0))[0]
        filtered_lining = np.take(raw_alpha_lining,filtered_idx,axis=0)
        self.alpha_atom_lining = filtered_lining
        filtered_xyz = np.take(raw_alpha_xyz,filtered_idx,axis=0)

        # cluster the remaining vertices to give index of belonging pockets
        zmat = linkage(filtered_xyz,method='average')
        cluster = fcluster(zmat,self.config.pocket_cluster_distance / 10,criterion='distance') # /10 turn A to nm
        self.alpha_atom_pocket_index = cluster -1 # cluster index start from 1

        # Reorganize into list of pockets
        self.pocket_alpha_atoms = [[] for _ in range(max(self.alpha_atom_pocket_index)+1)]
        for alpha_cluster_i, alpha_atom_idx in sorted(zip(self.alpha_atom_pocket_index, range(len(self.alpha_atom_pocket_index)))):
            self.pocket_alpha_atoms[alpha_cluster_i].append(alpha_atom_idx)

        for _ in self.pocket_alpha_atoms[1:]:
            self.top.add_residue(name='ASC',chain=self.top.chain(0))
        for i,pocket_index in enumerate(self.alpha_atom_pocket_index):
            self.top.add_atom('AAC',None,self.top.residue(self.alpha_atom_pocket_index[i]),i)

        self.traj.xyz = np.expand_dims(filtered_xyz,axis=0)

        # todo generate bac and acc and set as pocket residue atom name
        # todo write methods for retrieving these information

        filtered_lining_xyz = np.take(self.receptor_traj.xyz[0],self.alpha_atom_lining,axis=0)

        # calculate the polarity of alpha atoms
        total_score = np.array([getTetrahedronVolume(i) for i in filtered_lining_xyz])
        polar_ratio = np.average(np.take(self.atom_polarity,self.alpha_atom_lining),axis=1)
        polar_score = total_score * polar_ratio
        nonpolar_score = total_score - polar_score

        self.alpha_atom_score = np.stack([total_score, polar_score, nonpolar_score], axis=-1)*1000
        # also initialize contact score
        self.alpha_atom_contact = np.zeros(self.top.n_atoms,dtype=int)

    def build_beta_atoms(self):
        """
        cluster all alpha atoms into either beta atoms or pocket atoms (gamma)
        :param into: str, 'gamma' or 'beta'
        :return:
        """
        zmat = linkage(self.traj.xyz[0],method='average')
        cluster = fcluster(zmat,self.config.beta_cluster_cutoff / 10,criterion='distance')

    def calculate_contact(self,coordinates):
        """
        calculate which point in the alpha cluster is in contact with the given coordinates
        :param coordinates: Array
        :return:
        """
        np.put(self.alpha_atom_contact,
               checkContact(self.traj.xyz[0],coordinates,threshold=self.config.contact_threshold / 10)[0],1)
        self.alpha_atom_space_contact = np.transpose(np.transpose(self.alpha_atom_score) * self.alpha_atom_contact)

    def get_lining_atoms(self, index_list):
        """
        get a set of surface lining atoms in the given cluster of alpha atoms
        :type index_list: list
        """
        return np.unique(np.take(self.alpha_atom_lining, index_list).flatten())

    def get_lining_residues(self, index_list):
        # lining_atom_idx = np.take(self.alpha_atom_lining, index_list).flatten()

        lining_atom = [self.receptor_top.atom(i) for i in self.get_lining_atoms(index_list)]
        lining_residue_idx = [atom.residue.index for atom in lining_atom]

        return np.unique(lining_residue_idx)

    def get_cluster_centroid(self,cluster):
        xyz = np.take(self.traj.xyz[0],cluster,axis=0)
        return np.mean(xyz,axis=0)


    def screen_by_contact(self):

        contact_alpha_atoms = np.where(self.alpha_atom_contact > 0 )[0]
        contact_pocket = np.unique(np.take(self.alpha_atom_pocket_index,contact_alpha_atoms))
        #
        # contact_pocket_alpha_atom_index = [np.array(np.where(self.alpha_atom_pocket_index == i)) for i in contact_pocket]

        contact_pocket_alpha_atom_index = np.array(np.where(self.alpha_atom_pocket_index == contact_pocket[0]))
        for i in contact_pocket[1:]:
            contact_pocket_alpha_atom_index = np.append(contact_pocket_alpha_atom_index,np.array(np.where(self.alpha_atom_pocket_index == i)))
        self._update_alpha_atoms(contact_pocket_alpha_atom_index)

    def _update_alpha_atoms(self, new_alpha_index):
        self.traj.atom_slice(new_alpha_index, inplace=True)
        self.top = self.traj.top

        for pocket_idx,pocket in enumerate(self.top.residues):
            pocket.index = pocket_idx
        for atom_idx,atom in enumerate(self.top.atoms):
            atom.index = atom_idx

        self.alpha_atom_pocket_index = [atom.residue.index for atom in self.top.atoms]
        self.alpha_atom_lining = np.take(self.alpha_atom_lining, new_alpha_index, axis=0)
        self.alpha_atom_contact = np.take(self.alpha_atom_contact, new_alpha_index,axis=0)
        self.alpha_atom_score = np.take(self.alpha_atom_score, new_alpha_index, axis=0)

        self.pocket_alpha_atoms = [[] for _ in range(max(self.alpha_atom_pocket_index) + 1)]
        for alpha_cluster_i,alpha_atom_idx in sorted(
                zip(self.alpha_atom_pocket_index,range(len(self.alpha_atom_pocket_index)))):
            self.pocket_alpha_atoms[alpha_cluster_i].append(alpha_atom_idx)

        for pocket_idx,pocket in enumerate(self.top.residues):
            pocket.index = pocket_idx
        for atom_idx,atom in enumerate(self.top.atoms):
            atom.index = atom_idx


    def _build_pockets(self):
        pockets = []
        for p_residue in self.top.residues:
            alpha_atom_idx = self.pocket_alpha_atoms[p_residue.index]
            p_total_score = np.sum(np.take(self.alpha_atom_score, alpha_atom_idx, axis=0), axis=0)
            p_contact_score = np.sum(
                np.take(self.alpha_atom_score * np.expand_dims(self.alpha_atom_contact, axis=1),
                        alpha_atom_idx, axis=0), axis=0)
            lining_atoms = self.get_lining_atoms(alpha_atom_idx)
            lining_residues = self.get_lining_residues(alpha_atom_idx)
            pocket = AS_Pocket(contact_score=p_contact_score, total_score=p_total_score,parent=self, index=p_residue.index)
            pocket.lining_atoms = lining_atoms
            pocket.lining_residues = lining_residues
            pockets.append(pocket)
        self.pockets = pockets


class AS_Pocket:
    def __init__(self,contact_score,total_score,parent,index = 0):
        self.index = index
        self.parent = parent
        self.snapshot_idx = self.parent.snapshot_idx
        self._contact_score = contact_score
        self._total_score = total_score
        self._component_lookup = {'all':0, 'polar':1, 'nonpolar':2}


    def get_contact_score(self, component = 'all'):
        if component in self._component_lookup:
            return self._contact_score[self._component_lookup[component]]
        else:
            raise("Invalid component name, please use all, polar, or nonpolar, you used {}".format(component))

    def get_total_score(self,component = 'all'):
        if component in self._component_lookup:
            return self._total_score[self._component_lookup[component]]
        else:
            raise ("Invalid component name, please use all, polar, or nonpolar, you used {}".format(component))

    def __repr__(self):
        return "Pocket {} with total score of {}".format(self.index,self.get_total_score('all'))