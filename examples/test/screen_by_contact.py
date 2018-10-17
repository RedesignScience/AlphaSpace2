"""
Example:

    How to mark and select pocket points by contact

"""

import alphaspace
import sys
import mdtraj as md

if __name__ == '__main__':
    u = alphaspace._load(sys.argv[1])

    ligand_traj = md.load(sys.argv[2])

    u.calculateContact(ligand_traj.xyz[0])

    for snapshot in u.snapshots.values():

        for i, is_contact in enumerate(snapshot.pocket_contact):
            print(i, is_contact)
