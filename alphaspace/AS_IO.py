from __future__ import print_function, division

from multiprocessing import Pool

import mdtraj as md
import numpy as np
from alphaspace.AS_Snapshot import AS_Snapshot
from alphaspace.AS_Universe import AS_Universe
from alphaspace.AS_Vina import get_probe_score, pre_process_pdbqt
from alphaspace.AS_Funct import bin_cluster, group


def _get_path(f):
    """

    Parameters
    ----------
    f: list or str

    Returns
    -------

    """
    if isinstance(f, list):
        return f
    elif isinstance(f, str):
        with open(f, 'r') as handle:
            lines = handle.read().splitlines()
        return lines
    else:
        raise Exception


def _load_traj_with_ref(*args, top=None, ref=None):
    """

    Parameters
    ----------
    args: str
    top: str
    ref: str

    Returns
    -------

    """

    trajectory = []
    for arg in args:
        traj = md.load_pdb(arg) if arg.endswith('pdbqt') else md.load(arg, top)
        trajectory.append(traj)
    trajectory = md.join(trajectory)

    if ref is not None and ref.endswith('pdbqt'):
        with open(ref, 'r') as handle:
            lines = [line for line in handle.read().splitlines() if line.startswith("ATOM")]
        atom_numbers = []
        partial_charges = []
        adv_atom_types = []

        for pdb_line in lines:
            serial_number = int(pdb_line[6:11])
            atom_numbers.append(serial_number)
            name_with_spaces = pdb_line[12:16]
            alternate_location_indicator = pdb_line[16]
            residue_name_with_spaces = pdb_line[17:20]
            residue_number = int(pdb_line[22:26])
            insertion_code = pdb_line[26]
            x = float(pdb_line[30:38])
            y = float(pdb_line[38:46])
            z = float(pdb_line[46:54])
            partial_charge = float(pdb_line[70:76])
            partial_charges.append(partial_charge)
            adv_atom_type = pdb_line[77:79].strip()
            adv_atom_types.append(adv_atom_type)

        retaining_atom_indices = [atom.index for atom in trajectory.top.atoms if
                                  atom.serial in set(atom_numbers)]
        trajectory.atom_slice(atom_indices=retaining_atom_indices, inplace=True)

        trajectory.partial_charges = partial_charges
        trajectory.adv_atom_types = adv_atom_types


    else:

        print(
            'No reference structure selected, using all atoms in trajectory file and no vina score will be calculated'
        )
        trajectory.partial_charges = None
        trajectory.adv_atom_types = None

    return trajectory


def load(f, top=None, ref=None):
    file_path = _get_path(f)
    trajectory = _load_traj_with_ref(*file_path, top=top, ref=ref)
    return trajectory


def process(trajectory):
    snapshots = []
    for i in range(trajectory.n_frames):
        snapshot = trajectory[i]
        snapshot.adv_atom_types = trajectory.adv_atom_types
        snapshot.partial_charges = trajectory.partial_charges
        snapshots.append(snapshot)

    u = AS_Universe()

    # results  = []
    # for isnapshot in enumerate(snapshots):
    #     results.append(_process_snapshot(isnapshot))
    # #
    with Pool() as pool:
        results = pool.map(_process_snapshot, enumerate(snapshots))
    #     pool.join()
    #     pool.close()



    u.update(dict(results))

    return u


def _process_snapshot(isnapshot):
    print(isnapshot)
    i, snapshot = isnapshot
    ss = AS_Snapshot()

    print(snapshot)

    ss.tessellation(snapshot, snapshot_idx=0)
    ss._gen_beta()
    ss._gen_pocket()

    if snapshot.adv_atom_types is not None:
        prot_types, hp_type, acc_type, don_type = pre_process_pdbqt(snapshot)
        ss.beta_scores = get_probe_score(probe_coords=ss.beta_xyz * 10, prot_coord=snapshot.xyz[0] * 10,
                                         prot_types=prot_types,
                                         hp_type=hp_type,
                                         acc_type=acc_type, don_type=don_type)
    else:
        ss.beta_scores = np.zeros(len(ss.beta_xyz), dtype=np.float)
        print("No Vina Typing found, this is because you did not supply a pdbqt file. No Beta Score will be computed")
    return i, ss


# def _load_and_run_score(i_item):
#     i, item = i_item
#
#     no_pdbqt = True
#     if isinstance(item, str):
#         if len(item) >= 5 and item[-5:] == 'pdbqt':
#             traj = load_pdbqt(item)
#             no_pdbqt = False
#         else:
#             traj = md.load(item)
#     else:
#         traj = item
#
#     traj = traj.atom_slice(atom_indices=[atom.index for atom in traj.top.atoms if atom.element.atomic_number != 1])
#     ss = AS_Snapshot()
#     ss.tessellation(traj, 0)
#     ss._gen_beta()
#     ss._gen_pocket()
#
#     if no_pdbqt:
#         prot_types, hp_type, acc_type, don_type = gen_vina_type(ss.atom_names, ss.residue_names, ss.elements)
#     else:
#         prot_types, hp_type, acc_type, don_type = pre_process_pdbqt(traj)
#
#     ss.beta_scores = get_probe_score(probe_coords=ss.beta_xyz * 10, prot_coord=traj.xyz[0] * 10, prot_types=prot_types,
#                                      hp_type=hp_type,
#                                      acc_type=acc_type, don_type=don_type)
#
#     return i, ss


# def _load_and_run(i_item):
#     i, item = i_item
#
#     if isinstance(item, str):
#         traj = md.load(item)
#     else:
#         traj = item
#     traj = traj.atom_slice(atom_indices=[atom.index for atom in traj.top.atoms if atom.element.atomic_number != 1])
#     ss = AS_Snapshot()
#     ss.tessellation(traj, 0)
#     ss._gen_beta()
#     ss._gen_pocket()
#
#     ss.beta_scores = np.zeros(len(ss.beta_xyz), dtype=np.float)
#
#     return i, ss
#
# def _load(files, score=True):
#     """
#
#     Parameters
#     ----------
#     files: list
#
#     Returns
#     -------
#
#     """
#     if type(files) is not list:
#         files = [files]
#     u = AS_Universe()
#
#     if score:
#         with Pool() as pool:
#             results = pool.map(_load_and_run_score, enumerate(files), 1)
#     else:
#         with Pool() as pool:
#             results = pool.map(_load_and_run, enumerate(files), 1)
#     # print(results)
#     u.update(dict(results))
#
#     return u


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
