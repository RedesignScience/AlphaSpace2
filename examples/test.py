import alphaspace as al
import mdtraj
import os

# Load in receptor and ligand seperately using mdtraj. So you can load anything that is supported by mdtraj.
receptor = mdtraj.load("mdm2_p53/mdm2.pdb")
binder = mdtraj.load("mdm2_p53/p53.pdb")

# If you have the pdbqt file of the receptor, you can calculate the Vina Score. You can load the pdbqt file like this
al.annotateVinaAtomTypes(pdbqt="mdm2_p53/mdm2.pdbqt", receptor=receptor)

# Initialize a snapshot object, this will contain the receptor and the binder informations
ss = al.Snapshot()
ss.beta_cluster_dist = 1.6
ss.contact_cutoff = 1.6
ss.pocket_cluster_dist = 4.7
# Run the snapshot object by feeding it receptor and binder mdtraj objects.
ss.run(receptor=receptor, binder=binder)


with open('test.pdb','w') as handle:
	al.write_to_pdb(ss,handle=handle,)

