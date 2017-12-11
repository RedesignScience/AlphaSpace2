

from openeye.oechem import *
from openeye.oedocking import *
import sys



if __name__ == '__main__':
    
    imstr = oemolistream(sys.argv[1])
    proteinStructure = OEGraphMol()
    OEReadMolecule(imstr, proteinStructure)

    box = OEBox(-30	, 60, 32, -40, 40, 20)
    receptor = OEGraphMol()
    OEMakeReceptor(receptor, proteinStructure, box)