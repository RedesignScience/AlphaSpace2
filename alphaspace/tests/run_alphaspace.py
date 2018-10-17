import mdtraj
import alphaspace as al


def run(input_file, lig_file, output_file, config_path=None):
    pass


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", dest="input_file", required=True,
                        help="input coordinate or trajectory file", metavar="FILE", )
    parser.add_argument("-l", dest="ligand_file", required=False,
                        help="input for ligand file", metavar="FILE", )

    parser.add_argument("-o", dest="output_file", required=True,
                        help="output file, .pdb text or .as binary", metavar="FILE", )
    parser.add_argument("-c", dest="config_file", required=False,
                        help="configuration file", metavar="FILE")
    parser.add_argument("--verbose", help="increase output verbosity", required=False,
                        action="store_true")
    args = parser.parse_args()

    run(input_file=args.input_file, lig_file=args.ligand_file, output_file=args.output_file,
        config_path=args.config_file)
