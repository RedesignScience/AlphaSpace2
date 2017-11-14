import mdtraj

import alphaspace
from alphaspace.AS_IO import get_pdb_line


def run_alphaspace(pdb):
    universe = alphaspace.AS_Universe(receptor=mdtraj.load(pdb))

    universe.run_alphaspace()

    universe.screen_pockets()

    "ATOM      1  N   GLU     1      10.801 -12.147  -5.180  1.00  0.00"
    "ATOM       alpha_idx  A   AAC     Pocket_idx      xyz  polar  nonpolar polar_occupy_rate nonpolar_occupy_rate"
    "ATOM       pocket_center  A   ACC     Pocket_idx      xyz  polar  nonpolar polar_occupy_rate nonpolar_occupy_rate"
    for pocket in universe.pockets():
        for alpha in pocket.alphas:
            print(get_pdb_line(atomnum=alpha.idx, atomname='A', resname="AAC", resnum=pocket._idx, x=alpha.xyz[0],
                               y=alpha.xyz[0], z=alpha.xyz[0], occ=alpha.polar_score, bfactor=alpha.nonpolar_score))

        print(get_pdb_line(atomnum=alpha.idx, atomname='A', resname="AAC", resnum=pocket._idx, x=alpha.xyz[0],
                           y=alpha.xyz[0], z=alpha.xyz[0], occ=alpha.polar_score, bfactor=alpha.nonpolar_score))

    return


if __name__ == '__main__':

    import sys

    pdb = sys.argv[1]

    if len(sys.argv[1:]) == 2:
        output_dir = sys.argv[2]
    else:
        output_dir = '.'

    run_alphaspace(pdb)
