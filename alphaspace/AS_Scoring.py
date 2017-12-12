"""
This is the scoring module of AlphaSpace,
it aims to integrate multiple tool kits for easy analysis.
"""

from openeye.oechem import *
from openeye.oedocking import *
import numpy as np


def score_ligand(*args):


    receptor = OEGraphMol()

    coord =np.array([a[1] for a in args[0].GetCoords().items()])

    max_xyz = coord.max(axis=0) + 10
    min_xyz = coord.min(axis=0) - 10

    box = OEBox(max_xyz[0],max_xyz[1],max_xyz[2],min_xyz[0],min_xyz[1],min_xyz[2])

    print('box generated')

    OEMakeReceptor(receptor, args[0], box)


    print('box generated')

    # dock = OEDock(OEDockMethod_Chemgauss4,OESearchResolution_Default)
    # dock.Initialize(receptor)

    score = OEScore(OEScoreType_Chemgauss3)
    score.Initialize(receptor)


    print(score.ScoreLigand(args[1]))


    def PrintScore(dock, pose):
        print("Total pose score = %f" % dock.ScoreLigand(pose))
        print("Score components contributions to score:")
        for comp in dock.GetComponentNames():
            print("%12s: %6.2f" % (comp, dock.ScoreLigandComponent(pose, comp)))


    def PrintAtomScore(dock, pose, atom):
        print("")
        print("  Atom: %d score: %f" % (atom.GetIdx(), dock.ScoreAtom(atom, pose)))
        print("Score components contributions to atoms score: ")
        for comp in dock.GetComponentNames():
            print("%12s: %.2f" % (comp, dock.ScoreAtomComponent(atom, pose, comp)))

    # PrintScore(score,pose=args[1])


def gen_ligand(ligand: OEGraphMol,element):
    new_lig = OEGraphMol()


    coord = [a[1] for a in ligand.GetCoords().items()]

    for i in range(ligand.NumAtoms()):
        new_lig.NewAtom(element)

    new_lig.SetCoords(np.concatenate(coord))

    return new_lig






if __name__ == '__main__':
    import sys

    # ifs = oemolistream()
    # ofs = oemolostream()
    #
    # ifs.open(sys.argv[1])
    # prot_mol = next(ifs.GetOEGraphMols())
    #
    # ifs.open(sys.argv[2])
    # lig_mol  = next(ifs.GetOEGraphMols())
    #
    # # lig_mol = gen_ligand(lig_mol,6)
    #
    # print('structure loaded')
    # score_ligand(prot_mol,lig_mol)
    #
    #

    imstr = oemolistream(sys.argv[1])
    proteinStructure = OEGraphMol()
    OEReadMolecule(imstr, proteinStructure)
    box = OEBox(-35.959+5.77804, 52.7634 + 7 , 26.3825 + 7 ,-35.959-5.77804,  52.7634 - 7 , 26.3825 - 7 )

    # coord = np.array([a[1] for a in proteinStructure.GetCoords().items()])
    #
    # max_xyz = coord.max(axis=0) + 10
    # min_xyz = coord.min(axis=0)
    #
    # box = OEBox(max_xyz[0], max_xyz[1], max_xyz[2], min_xyz[0], min_xyz[1], min_xyz[2])

    receptor = OEGraphMol()
    OEMakeReceptor(receptor, proteinStructure, box)