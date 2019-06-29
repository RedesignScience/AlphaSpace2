"""
To provide compatibility to original AlphaSpace, write out files to a designated folder, and create visualization file.
"""

from .Cluster import _Alpha, _Beta, _Pocket

import os
import shutil

def _format_83(f):
    """Format a single float into a string of width 8, with ideally 3 decimal
    places of precision. If the number is a little too large, we can
    gracefully degrade the precision by lopping off some of the decimal
    places. If it's much too large, we throw a ValueError"""
    if -999.999 < f < 9999.999:
        return '%8.3f' % f
    if -9999999 < f < 99999999:
        return ('%8.3f' % f)[:8]
    raise ValueError('coordinate "%s" could not be represented '
                     'in a width-8 field' % f)


def write_chimera_scripts(folder):
    from shutil import copy

    import os

    _ROOT = os.path.abspath(os.path.dirname(__file__))

    def get_data(path):
        return os.path.join(_ROOT, 'data', path)
    if not os.path.isdir(folder):
        os.mkdir(folder)

    copy(get_data('AS_Chimera_beta.py'), folder)
    copy(get_data('AS_Chimera_alpha.py'), folder)
    copy(get_data('colors_table.txt'), folder)
    copy(get_data('colors_chimera.txt'), folder)


def write_snapshot(folder_path, snapshot, receptor=None, binder=None, chimera_scripts=True, contact_only=True):


    if os.path.isdir(os.path.join(folder_path, 'pockets')):
        shutil.rmtree(os.path.join(folder_path, 'pockets'))
    os.makedirs(os.path.join(folder_path, 'pockets'))


    if chimera_scripts:
        write_chimera_scripts(folder_path)

    if receptor or binder:
        if not os.path.isdir(os.path.join(folder_path, 'pdb_out')):
            os.makedirs(os.path.join(folder_path, 'pdb_out'))
    if receptor:
        receptor.save(os.path.join(folder_path, 'pdb_out', 'prot.pdb'))
    if binder:
        binder.save(os.path.join(folder_path, 'pdb_out', 'lig.pdb'))

    pocket_index = 0

    if contact_only:

        pockets = sorted([p for p in snapshot.pockets if p.isContact], key=lambda p: p.space, reverse=True)
    else:
        pockets = sorted([p for p in snapshot.pockets], key=lambda p: p.space, reverse=True)

    for pocket in pockets:
        pocket_index += 1
        if receptor:
            lining_atoms = receptor.atom_slice(pocket.lining_atoms_idx)
            lining_atoms.save(os.path.join(folder_path, 'pockets', '{}_alpha.pdb'.format(pocket_index)))
            lining_atoms.save(os.path.join(folder_path, 'pockets', '{}_beta.pdb'.format(pocket_index)))

        with open(os.path.join(folder_path, 'pockets', '{}_beta.pdb'.format(pocket_index)), 'a') as handle:
            for beta in pocket.betas:
                handle.write(gen_pdb_line(atomIndex=beta.index,
                                          atomName='BAO' if beta.isContact else 'BAU',
                                          resName='BAC',
                                          resIndex=pocket_index,
                                          chainName=" ",
                                          bfactor=beta.score,
                                          element=beta.best_probe_type,
                                          xyz=beta.centroid,
                                          occupancy=" ")
                             )
            handle.write(gen_pdb_line(atomIndex=pocket_index,
                                      atomName='BCC',
                                      resName='BCC',
                                      resIndex=pocket_index,
                                      chainName=" ",
                                      bfactor=pocket.score,
                                      element='C',
                                      xyz=pocket.centroid))

        with open(os.path.join(folder_path, 'pockets', '{}_alpha.pdb'.format(pocket_index)), 'a') as handle:
            for alpha in pocket.alphas:
                line = gen_pdb_line(atomIndex=alpha.index,
                                    atomName='AAO' if alpha.isContact else 'AAU',
                                    resName='AAC',
                                    resIndex=pocket_index,
                                    chainName=" ",
                                    bfactor=alpha.space,
                                    occupancy=" ",
                                    element='C',
                                    xyz=alpha.centroid,
                                    )
                handle.write(line)
            handle.write(gen_pdb_line(atomIndex=pocket_index,
                                      atomName='ACC',
                                      resName='ACC',
                                      resIndex=pocket_index,
                                      chainName=" ",
                                      bfactor=pocket.score,
                                      element='C',
                                      xyz=pocket.centroid))



def gen_pdb_line(atomIndex, atomName, resName, resIndex, chainName, bfactor, element, xyz, occupancy=" "):
    line = "ATOM  %5d %-4s %3s %1s%4d    %s%s%s  1.00 %5.2f      %-4s%-2s  \n" % (
        atomIndex % 100000, atomName, resName, chainName,
        resIndex % 10000, _format_83(xyz[0]),
        _format_83(xyz[1]), _format_83(xyz[2]),
        bfactor, occupancy, element)
    return line
