import mdtraj

import alphaspace


def run_alphaspace(input_file, output_file, config=None):
    universe = alphaspace.AS_Universe(receptor=mdtraj.load(input_file))


    universe.run_alphaspace()

    universe.screen_pockets()

    if output_file.endswith('.as'):
        """save as .as binary"""
        alphaspace.save_universe(universe, output_file)

    elif output_file.endswith('.pdb'):
        """save as pdb readable file"""
        out_lines = []
        for pocket in universe.pockets():
            for alpha in pocket.alphas:
                out_lines.append(
                    alphaspace.get_pdb_line(atomnum=alpha.idx, atomname='A', resname="AAC", resnum=pocket._idx,
                                            x=alpha.xyz[0],
                                            y=alpha.xyz[0], z=alpha.xyz[0], occ=alpha.polar_score,
                                            bfactor=alpha.nonpolar_score))
        with open(output_file, 'r') as handle:
            handle.writelines(out_lines)


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", dest="input_file", required=True,
                        help="input coordinate or trajectory file", metavar="FILE", )
    parser.add_argument("-o", dest="output_file", required=True,
                        help="output file, .pdb text or .as binary", metavar="FILE", )
    parser.add_argument("-c", dest="config_file", required=False,
                        help="configuration file", metavar="FILE")
    parser.add_argument("--verbose", help="increase output verbosity", required=False,
                        action="store_true")
    args = parser.parse_args()

    run_alphaspace(input_file=args.input_file, output_file=args.output_file, config=args.config_file)
