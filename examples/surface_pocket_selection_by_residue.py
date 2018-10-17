"""
This file shows you how you can select pockets by screening them against a list of anchoring residues
"""

import alphaspace
import mdtraj
import sys

receptor_path, binder_path = sys.argv[1], sys.argv[2]

universe = alphaspace.Trajectory()
universe.set_receptor(structure=mdtraj.load(receptor_path), keepH=True)
universe.set_binder(structure=mdtraj.load(binder_path))

universe.run_alphaspace()

for i, residue in enumerate(universe.receptor.top.residues):
    print(i, residue)

for residue in universe.receptor.top.residues:
    print(residue)

# let's say you only want pockets that have lining atoms from residue 100,101,102, you can choose to deactivate all
# that does not satisfy this criteria.

for pocket in universe.pockets(snapshot_idx=0, active_only=False):
    if not set(pocket.lining_residues_idx).intersection([100, 101, 102]):
        pocket.deactivate()
