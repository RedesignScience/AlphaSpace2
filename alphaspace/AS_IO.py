from __future__ import print_function, division

import pickle
from multiprocessing import Pool

import numpy as np
import mdtraj as md
from alphaspace.AS_Snapshot import AS_Snapshot
from alphaspace.AS_Universe import AS_Universe
from alphaspace.AS_Vina import gen_vina_type, get_probe_score, pre_process_pdbqt


def _load(files):
    """

    Parameters
    ----------
    files: list

    Returns
    -------

    """
    if type(files) is not list:
        files = [files]
    u = AS_Universe()

    with Pool() as pool:
        results = pool.map(_load_and_run, enumerate(files), 1)

    print(results)
    # print(results)
    u.update(dict(results))

    return u


def _load_and_run(i_item):
    i, item = i_item

    if isinstance(item, str):
        if len(item) >= 5 and item[-5:] == 'pdbqt':
            traj = load_pdbqt(item)
            prot_types, hp_type, acc_type, don_type = pre_process_pdbqt(traj)
        else:
            traj = md.load(item)
    else:
        traj = item

    traj = traj.atom_slice(atom_indices=[atom.index for atom in traj.top.atoms if atom.element.atomic_number != 1])
    ss = AS_Snapshot()
    ss.tessellation(traj, 0)
    ss._gen_beta()
    ss._gen_pocket()

    # prot_types, hp_type, acc_type, don_type = gen_vina_type(ss.atom_names,
    #                                                                              ss.residue_names,
    #                                                                              ss.elements)

    ss.beta_scores = get_probe_score(probe_coords=ss.beta_xyz * 10, prot_coord=traj.xyz[0] * 10, prot_types=prot_types,
                                     hp_type=hp_type,
                                     acc_type=acc_type, don_type=don_type)
    return i, ss


def load_from_file(path_file):
    with open(path_file, 'r') as handle:
        lines = handle.read().splitlines()
    return _load([line for line in lines if (not line.startswith('#')) and len(line) > 1])


def load_pdbqt(filename, stride=None, atom_indices=None, frame=None, no_boxchk=False, standard_names=False):
    traj = md.load_pdb(filename, stride, atom_indices, frame,
                       no_boxchk, standard_names)

    traj.atom_slice([atom.index for atom in traj.top.atoms if atom.element.atomic_number != 1], inplace=True)

    charges = []
    pdbqt_atom_name = []
    with open(filename, 'r') as f:
        for line in f.read().splitlines():
            if line[:6] in {'ATOM  ', 'HETATM'}:
                row = line.split()
                pdbqt_name = str(row[-1]).strip()
                charge = float(str(row[-2]).strip())
                if pdbqt_name not in {'H', 'HD', 'HS'}:
                    charges.append(charge)
                    pdbqt_atom_name.append(pdbqt_name)

    for atom in traj.top.atoms:
        atom.charge = charges[atom.index]
        atom.pdbqt_name = pdbqt_atom_name[atom.index]
    return traj


def get_pdb_line(atom='ATOM', atomnum=0, atomname=' ', resname=' ', chain_idx=' ', resnum=0, x=0.0, y=0.0, z=0.0,
                 occ=0.0,
                 bfactor=0.0, element=' '):
    j = [str(atom), str(int(atomnum)), str(atomname), str(resname), str(chain_idx), str(int(resnum)), float(x),
         float(y), float(z), float(occ), float(bfactor), str(element)]
    j[0] = j[0].ljust(6)  # atom#6s
    j[1] = j[1].rjust(5)  # aomnum#5d
    j[2] = j[2].center(4)  # atomname$#4s
    j[3] = j[3].ljust(3)  # resname#1s
    j[4] = j[4].rjust(1)  # Astring
    j[5] = j[5].rjust(4)  # resnum
    j[6] = str('%8.3f' % j[6]).rjust(8)  # x
    j[7] = str('%8.3f' % j[7]).rjust(8)  # y
    j[8] = str('%8.3f' % j[8]).rjust(8)  # z\
    j[9] = str('%6.2f' % (j[9])).rjust(6)  # occ
    j[10] = str('%6.2f' % j[10]).ljust(6)  # temp
    j[11] = j[11].rjust(12)  # elname
    line = ("{}{} {} {} {}{}    {}{}{}{}{}{}".format(
        j[0], j[1], j[2], j[3], j[4], j[5], j[6], j[7], j[8], j[9], j[10], j[11]))

    return line


def save_universe(universe, file_path: str) -> bool:
    try:
        with open(file_path, 'wb') as handle:
            pickle.dump(universe, handle)
        return True
    except:
        print('Cannot save to {}'.format(file_path))
        return False


def load_universe(file_path: str):
    try:
        with open(file_path, 'rb') as handle:
            return pickle.load(handle)
    except:
        print('Cannot load from {}'.format(file_path))
