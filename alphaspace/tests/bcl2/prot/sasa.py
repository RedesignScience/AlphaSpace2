import numpy as np
import mdtraj as md
import sys

pdb = sys.argv[1]
out = sys.argv[2]

trajectory = md.load(pdb)
sasa = md.shrake_rupley(trajectory,n_sphere_points=2000,)



with open(out,'w') as handle:
	for atom in trajectory.topology.atoms:
		# print(out)
		handle.write("{} {} {} {}\n".format(atom.index,atom.name,atom.residue.index, sasa[0][atom.index]*100))
