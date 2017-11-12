import mdtraj

import alphaspace


def run_alphaspace(pdb):
    universe = alphaspace.AS_Universe(receptor=mdtraj.load(pdb))

    universe.run_alphaspace()

    universe.screen_pockets()

    for pocket in universe.pockets():
        print(len(pocket.alpha_idx), pocket.polar_score(), pocket.nonpolar_score())

    return


if __name__ == '__main__':

    import sys

    pdb = sys.argv[1]

    if len(sys.argv[1:]) == 2:
        output_dir = sys.argv[2]
    else:
        output_dir = '.'

    run_alphaspace(pdb)
