import sys, os
import numpy as np
import logging
from scipy import spatial
from scipy.spatial import cKDTree
from scipy.spatial.distance import cdist

polarTypes = np.array(['OA', 'O', 'N', 'NA', 'S', 'SA'])
nonPolarTypes = np.array(['C', 'A', 'Cl', 'CL', 'Br', 'BR', 'F', 'I'])

autodockVinaAtomTypes = {'H': 0.00, 'HD': 0.00, 'C': 1.90, 'A': 1.90,
                         'N': 1.80, 'NA': 1.80, 'OA': 1.70, 'O':1.70, 'S': 2.00,
                         'SA': 2.00, 'P': 2.10, 'F': 1.50, 'Cl': 1.80, 'CL': 1.80,
                         'Br': 2.00, 'BR': 2.00, 'I': 2.20, 'Mn': 1.20, 'MN': 1.20,
                         'Mg': 1.20, 'MG': 1.20, 'Zn': 1.20, 'ZN': 1.20, 'Ca': 1.20,
                         'CA': 1.20, 'Fe': 1.20, 'FE': 1.20, 'Ni': 1.20, 'NI': 1.20,
                         'Na': 1.20, 'Co': 1.20, 'CO': 1.20, 'Cu': 1.20, 'CU': 1.20,
                         'Cd': 1.20, 'CD': 1.20, 'K': 1.20, 'Se': 1.75, 'SE': 1.75}
                         
#autodockVinaTerms = np.array([-0.035579, -0.005156, 0.840245, -0.035069, -0.587439])
Lin_F9_Terms = np.array([-0.032879, -0.000896, -0.006915, 0.244318, -0.007459, -0.170027])
probElements = ['C', 'Br', 'F', 'Cl', 'I', 'OA', 'O', 'S', 'SA', 'N', 'NA', 'P']
cofactor_match_dict = {'C': 'C', 'N': 'N', 'P': 'P', 'O': 'OA', 'S': 'SA', 
                       'F': 'F', 'Cl': 'Cl', 'Br': 'Br','I': 'I'}

def _get_typing_dicts(hp_types_dat_path, autodock_atom_type_dat_path):
    hp_lines = np.loadtxt(hp_types_dat_path, dtype=str, delimiter=',')
    hp_types_dict = {str(i): {} for i in hp_lines[:, 0]}
    for i1, i2, i3, i4, i5 in hp_lines:
        hp_types_dict[i1][i2] = (i3, i4, i5)

    autodock_lines = np.loadtxt(autodock_atom_type_dat_path, dtype=str, delimiter=',')
    autodock_types_dict = {str(i).strip(): float(j) for i,j in autodock_lines}

    return hp_types_dict, autodock_types_dict


def _pre_process_pdb(pdb_lines):
    this_dir, this_filename = os.path.split(__file__)

    types_dict, autodock_types_dict = \
        _get_typing_dicts(os.path.join(this_dir, "LinF9data", "hp_types_dict.dat"),
                          os.path.join(this_dir, "LinF9data", "Vina_atom_radii.dat"))
                          
    polar_atoms = np.array(['O','OA','OS','N','NS','NA','S','SA','P'])
    ali_atoms = np.array(['C', 'A'])

    prot_types = []
    hp_type = []
    don_type = []
    acc_type = []
    prot_coord = []
    ### process heavy atom lines in pdb
    clean_pdb_lines = [t for t in pdb_lines if t[0:6] in ['ATOM  ', 'HETATM'] and t[76:].strip() not in ['H', 'HD', 'HS']]

    for ind, t in enumerate(clean_pdb_lines):
        resname = t[17:20].strip()
        atom_type = t[12:16].strip()
        if atom_type == 'OXT':  ### c terminal O is considered as OA
            atom_type = 'O'
        element_type = t[76:].strip()
        x_coords = float(t[30:38].strip())
        y_coords = float(t[38:46].strip())
        z_coords = float(t[46:54].strip())
        prot_coord.append([x_coords, y_coords, z_coords])
        ### prot_types is a list of autodock atom types for heavy atoms in pdb

        if resname in types_dict:
            if atom_type in types_dict[resname]:
                prot_types.append(types_dict[resname][atom_type][0])
        else:
            prot_types.append(cofactor_match_dict[element_type])
            
        if resname in types_dict and atom_type in types_dict[resname]:
            if prot_types[ind] in ali_atoms:
                hp_type.append(types_dict[resname][atom_type][1])
                don_type.append('XXX')
                acc_type.append('XXX')
            elif prot_types[ind] in polar_atoms:
                don_type.append(types_dict[resname][atom_type][2])
                acc_type.append(types_dict[resname][atom_type][1])
                hp_type.append('XXX')
            else:
                print(t)
                hp_type.append('UNK')
                don_type.append('UNK')
                acc_type.append('UNK')
        else:
            print(t)
            hp_type.append('UNK')
            don_type.append('UNK')
    prot_coord = np.array(prot_coord)
    prot_types = np.array(prot_types)
    hp_type = np.array(hp_type)
    acc_type = np.array(acc_type)
    don_type = np.array(don_type)

    return prot_types, hp_type, acc_type, don_type


def _pre_process_pdbqt(traj, truncation_length=0):
    def _assign_hp(prot_coord, prot_types, hp_type, tree=None):
        """
        """

        tree = cKDTree(prot_coord) if tree is None else tree

        for ix in np.where(hp_type == 'UNK')[0]:
            indx = np.array(tree.query_ball_point(prot_coord[ix], 2.0))
            if prot_types[ix] not in polarTypes:
                if np.any(np.in1d(prot_types[indx[indx != ix]], polarTypes)):
                    hp_type[ix] = 'NNP'
                else:
                    hp_type[ix] = 'NP'
            else:
                hp_type[ix] = 'XXX'
        return hp_type

    def _assign_acc(prot_types, acc_type):
        """
        """
        for ix in np.where(acc_type == 'UNK')[0]:
            if prot_types[ix] in ['OA', 'OS', 'SA', 'S']:
                acc_type[ix] = 'P'
            else:
                acc_type[ix] = 'XXX'
        return acc_type

    def _assign_don(prot_coord, prot_types, don_type, tree=None):
        tree = cKDTree(prot_coord) if tree is None else tree
        for ix in np.where(don_type == 'UNK')[0]:
            indx = np.array(tree.query_ball_point(prot_coord[ix], 2.0))
            if prot_types[ix] in ['N', 'NS', 'NA']:
                if len(indx[indx != ix]) > 2:
                    don_type[ix] = 'NPP'
                else:
                    don_type[ix] = 'P'
            else:
                don_type[ix] = 'XXX'
        return don_type
    
    this_dir, this_filename = os.path.split(__file__)
    #types_dict = autodockVinaAtomTypes  why????
    # change happens
    types_dict, autodock_types_dict = \
                _get_typing_dicts(os.path.join(this_dir, "LinF9data", "hp_types_dict.dat"),
                                  os.path.join(this_dir, "LinF9data", "Vina_atom_radii.dat"))
    ali_atoms = {'C', 'A'}

    hp_type = []
    don_type = []
    acc_type = []
    prot_coord = traj.xyz[0] * 10
    prot_tree = spatial.cKDTree(prot_coord)
    prot_types = []
    
    for i, atom in enumerate(traj.top.atoms):
        prot_types.append(traj.adv_atom_types[i])
        if atom.residue.name in types_dict.keys():
            if atom.name in types_dict[atom.residue.name].keys():
                ## if atom.pdbqt_name in ali_atom:  why???
                if types_dict[atom.residue.name][atom.name][0] in ali_atoms: ## change
                    hp_type.append(types_dict[atom.residue.name][atom.name][1])
                    don_type.append('XXX')
                    acc_type.append('XXX')
                elif types_dict[atom.residue.name][atom.name][0] in polarTypes: ## change
                    don_type.append(types_dict[atom.residue.name][atom.name][2])
                    acc_type.append(types_dict[atom.residue.name][atom.name][1])
                    hp_type.append('XXX')
                else:
                    hp_type.append('UNK')
                    don_type.append('UNK')
                    acc_type.append('UNK')
            else:
                hp_type.append('UNK')
                don_type.append('UNK')
                acc_type.append('UNK')
        else:
            hp_type.append('UNK')
            don_type.append('UNK')
            acc_type.append('UNK')

    prot_types = np.array(prot_types)
    hp_type = np.array(hp_type)
    acc_type = np.array(acc_type)
    don_type = np.array(don_type)
    
    if np.any(hp_type == 'UNK'):
        _assign_hp(prot_coord, prot_types, hp_type, prot_tree)
        _assign_don(prot_coord, prot_types, don_type, prot_tree)
        _assign_acc(prot_types, acc_type)

    if truncation_length != 0:
        return prot_types[:truncation_length], hp_type[:truncation_length], \
               acc_type[:truncation_length], don_type[:truncation_length], metal_type[:truncation_length]
    else:
        return prot_types, hp_type, acc_type, don_type


def _gen_vina_type(atom_names, residue_names, elements):
    """

    Parameters
    ----------
    pdb_lines

    Returns
    -------
    prot_types, hp_type, acc_type, don_type
    """

    this_dir, this_filename = os.path.split(__file__)

    types_dict, autodock_types_dict = _get_typing_dicts(
        os.path.join(this_dir, "LinF9data", "hp_types_dict.dat"),
        os.path.join(this_dir, "LinF9data", "Vina_atom_radii.dat"))

    polar_atoms = np.array(['OA', 'O', 'N', 'NA', 'S', 'SA'])
    ali_atoms = np.array(['C', 'A'])

    prot_types = []
    hp_type = []
    don_type = []
    acc_type = []
    ### process heavy atom lines in pdb

    for idx, atom_name, resname, element in zip(range(len(list(atom_names))), atom_names, residue_names, elements):
        if atom_name == 'OXT':  ### c terminal O is considered as OA
            atom_name = 'O'

        if resname in types_dict:
            if atom_name in types_dict[resname]:
                prot_types.append(types_dict[resname][atom_name][0])
        elif element in cofactor_match_dict:
            prot_types.append(cofactor_match_dict[element])
        else:
            raise Exception("{} if not a supported in vina scoring element".format(element))

        ### hp_type is acceptor type for aliphatic atoms
        #### don_type and acc_type is donor and acceptor type for polar atoms
        if resname in types_dict:
            if atom_name in types_dict[resname]:
                if prot_types[idx] in ali_atoms:
                    hp_type.append(types_dict[resname][atom_name][1])
                    don_type.append('XXX')
                    acc_type.append('XXX')
                elif prot_types[idx] in polar_atoms:
                    don_type.append(types_dict[resname][atom_name][2])
                    acc_type.append(types_dict[resname][atom_name][1])
                    hp_type.append('XXX')
            else:
                hp_type.append('UNK')
                don_type.append('UNK')
                acc_type.append('UNK')                
        else:
            hp_type.append('UNK')
            don_type.append('UNK')
            acc_type.append('UNK')

    prot_types = np.array(prot_types)
    hp_type = np.array(hp_type)
    acc_type = np.array(acc_type)
    don_type = np.array(don_type)
    return prot_types, hp_type, acc_type, don_type


def _get_probe_score(prot_coord, prot_types, hp_type, don_type, acc_type, probe_coords):
    """

    Examples
    --------

    Use beta atom as probe points and estimate the ligandibility
    Use the given protein atom type to calculate probe score at given location.

    Parameters
    ----------

    prot_coord : np.ndarray
        shape = (n,3)
    prot_types : np.ndarrayNPP
        shape = (n)
    hp_type : np.ndarray
        shape = (n)
    don_type : np.ndarray
        shape = (n)
    acc_type : np.ndarray
        shape = (n)
    probe_coords : np.ndarray
        shape = (m,3)

    Returns
    -------

    prb_dict : dict
        probe score dictionary


    """

    def _hydrophobic_interp(r):
        """
        """
        if r < 0.5:
            x = 1.0
        elif r > 1.5:
            x = 0.0
        else:
            x = 1.5 - r
        return x

    def _hbond_interp(r):
        if r < -0.7:
            x = 1.0
        elif r >= 0:
            x = 0.0
        else:
            x = -r / 0.7
        return x

    probe_prot_dist = cdist(probe_coords, prot_coord)

    probe_scores = np.zeros((probe_coords.shape[0], len(probElements)), dtype=np.float32)

    for probe_idx in range(probe_coords.shape[0]):

        dist_bool = probe_prot_dist[probe_idx] <= 12.0
        temp_dist = probe_prot_dist[probe_idx][dist_bool]
        
        temp_type = prot_types[dist_bool]
        NP_type = (hp_type[dist_bool] == 'NP')
        Pdon_type = (don_type[dist_bool] == 'P')
        Pacc_type = (acc_type[dist_bool] == 'P')
        dist_radii = np.array([autodockVinaAtomTypes[x] for x in temp_type])

        for element_idx, probe_element in enumerate(probElements):
            probe_dist = autodockVinaAtomTypes[probe_element]
            proc_dist = temp_dist - dist_radii - probe_dist
            g1 = np.sum(np.exp(-(proc_dist / 0.5) ** 2))
            #g2 = np.sum(np.exp(-((proc_dist - 3.0) / 2.0) ** 2))
            g61 = np.sum(np.exp(-((proc_dist - 2.5) / 0.5) ** 2))
            g62 = np.sum(np.exp(-((proc_dist - 5.0) / 0.5) ** 2))
            rep = np.sum([dd ** 2 if dd < 0.0 else 0.0 for dd in proc_dist])
            
            if probe_element in {'C', 'Br', 'Cl', 'F', 'I', 'S'}:
                h1 = np.sum([_hydrophobic_interp(dd) for dd in proc_dist[NP_type]])
                h2 = 0.0
            elif probe_element in {'OA', 'NA', 'SA'}:
                h1 = 0.0
                h2 = np.sum([_hbond_interp(dd) for dd in proc_dist[Pdon_type]])
            elif probe_element in {'N', 'O', 'P'}:
                h1 = 0.0
                h2 = np.sum([_hbond_interp(dd) for dd in proc_dist[Pacc_type]])
            else:
                raise Exception()

            probe_scores[probe_idx][element_idx] = np.sum(np.array([g1, g61, g62, rep, h1, h2]) * Lin_F9_Terms)

            # print(probe_scores[probe_idx][element_idx], g1,g2,rep,h1,h2)

    return probe_scores


def annotateVinaAtomTypes(receptor, pdbqt):
    """

    Parameters
    ----------
    receptor
    args: str
    top: str
    ref: str

    Returns
    -------

    """

    with open(pdbqt, 'r') as handle:
        lines = [line for line in handle.read().splitlines() if line.startswith(("ATOM","HETATM"))]
    atom_numbers = []
    partial_charges = []
    adv_atom_types = []

    for line in lines:
        serial_number = int(line[6:11])
        atom_numbers.append(serial_number)
        name_with_spaces = line[12:16]
        alternate_location_indicator = line[16]
        residue_name_with_spaces = line[17:20]
        residue_number = int(line[22:26])
        insertion_code = line[26]
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        partial_charge = float(line[70:76])
        partial_charges.append(partial_charge)
        adv_atom_type = line[77:79].strip()
        adv_atom_types.append(adv_atom_type)
        

    if receptor.top.n_atoms < len(adv_atom_types):
        print("Redundant Atom types in pdbqt file found, trimming last {} entries".format(
            {receptor.top.n_atoms - len(adv_atom_types)}))
        adv_atom_types = adv_atom_types[:receptor.top.n_atoms]
        partial_charges = partial_charges[:receptor.top.n_atoms]
    receptor.partial_charges = partial_charges
    receptor.adv_atom_types = adv_atom_types
