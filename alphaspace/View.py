"""
To provide compatibility to original AlphaSpace, write out files to a designated folder, and create visualization file.
"""

from .Cluster import _Alpha, _Beta, _Pocket

import os


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
    scripts_path = os.path.join(os.path.abspath(os.path.join(os.path.realpath(__file__), os.path.pardir)), 'data')

    if not os.path.isdir(folder):
        os.mkdir(folder)

    copy(os.path.join(scripts_path, 'AS_Chimera_beta.py'), folder)
    copy(os.path.join(scripts_path, 'AS_Chimera_alpha.py'), folder)
    copy(os.path.join(scripts_path, 'colors_table.txt'), folder)
    copy(os.path.join(scripts_path, 'colors_chimera.txt'), folder)


def write_snapshot(folder_path, snapshot, receptor=None):
    if not os.path.isdir(os.path.join(folder_path, 'pockets')):
        os.makedirs(os.path.join(folder_path, 'pockets'))

    for pocket in snapshot.pockets:
        if pocket.isContact:
            if receptor:
                lining_atoms = receptor.atom_slice(pocket.lining_atoms_idx)
                lining_atoms.save((folder_path, 'pockets', '{}_alpha.pdb'.format(pocket.index)))
                lining_atoms.save((folder_path, 'pockets', '{}_beta.pdb'.format(pocket.index)))

            with open(os.path.join(folder_path, 'pockets', '{}_beta.pdb'.format(pocket.index)), 'a') as handle:
                for beta in pocket.betas:
                    handle.write(gen_pdb_line(atomIndex=beta.index,
                                              atomName='BAO' if beta.isContact else 'BAU',
                                              resName='BAC',
                                              resIndex=pocket.index,
                                              chainName=" ",
                                              bfactor=beta.score,
                                              element='C',
                                              xyz=beta.centroid,
                                              occupancy=" ")
                                 )
                handle.write(gen_pdb_line(atomIndex=pocket.index,
                                          atomName='BCC',
                                          resName='BCC',
                                          resIndex=pocket.index,
                                          chainName=" ",
                                          bfactor=pocket.score,
                                          element='C',
                                          xyz=pocket.centroid))

            with open(os.path.join(folder_path, 'pockets', '{}_alpha.pdb'.format(pocket.index)), 'a') as handle:
                for alpha in pocket.alphas:
                    line = gen_pdb_line(atomIndex=alpha.index,
                                        atomName='AAO' if alpha.isContact else 'AAU',
                                        resName='AAC',
                                        resIndex=pocket.index,
                                        chainName=" ",
                                        bfactor=alpha.space,
                                        occupancy=" ",
                                        element='C',
                                        xyz=alpha.centroid,
                                        )
                    handle.write(line)
                handle.write(gen_pdb_line(atomIndex=pocket.index,
                                          atomName='ACC',
                                          resName='ACC',
                                          resIndex=pocket.index,
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
