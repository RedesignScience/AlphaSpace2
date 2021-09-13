import mdtraj as md
import alphaspace2 as al
from alphaspace2.Community import genCommunityPocket, extendedCommunityAnalysis

protein = md.load('your_protein.pdb')
protein = protein.atom_slice(protein.top.select("type != H")) # BE SURE TO STRIP HYDROGENS BEFORE RUNNING ALFASPACE!

snapshot = al.Snapshot()
snapshot.run(protein)

# Run community analysis saving alfa and beta atoms to separate PDB files prefixed with "protein." The files generated are protein_community_alfas.pdb, 
# protein_community_betas.pdb, and protein.npy (nonpolar space, [Å³], volume, [Å³], occluded ASA, [Å²] for community 1, nonpolar space, [Å³], volume, [Å³], occluded 
# ASA, [Å²] for community 2, etc. The communities are saved as separate residues (by residue index; core pocket atoms and auxilliary pocket atoms differ by name -
# A(B)CP or A(B)AP) with nonpolar spaces as b factor columns for easy visualization. community_properties is a tuple of thee elements - the first is a list of 
# community properties (the same ones saved in protein.npy), the second is a dictionary of the communities and their constituent pockets, the third is a list of 
# all pockets in the snapshot. 
   
community_properties = extendedCommunityAnalysis(protein = protein, snapshot = snapshot, atoms_to_writeout = 'both', save_npy = True, output_prefix = 'protein')
