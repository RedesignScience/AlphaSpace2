import numpy as np
import scipy.cluster.hierarchy as hier

from ..PDB.PDBParser import PDBParser


def fill_p_verts(p_verts, cluster):
    """
    Function: fill an empty list of lists (length equal to the number of clusters) with vertice indices according to the clustering (becomes centroids in super clustering)
    Parameters
    p_verts: empty list of lists
    cluster: cluster indices
    """

    for i, clust_idx in enumerate(cluster):
        p_idx = clust_idx - 1
        p_verts[p_idx].append(i)
    return


def cluster_coord_list(list_of_coord, cluster_method='average', cluster_distance=4.0):
    zmat = hier.linkage(list_of_coord, method=cluster_method)
    cluster = hier.fcluster(zmat, cluster_distance, criterion='distance')
    clusters = [[] for i in range(max(cluster))]
    fill_p_verts(clusters, cluster)
    frag_centroid = []
    for item in clusters:
        frag_centroid.append(np.average([list_of_coord[index] for index in item], axis=0))
    return frag_centroid


def gmc_coord_list(list_of_coord):
    return np.average(list_of_coord,axis=0)


def get_lfc(in_file,mode,out_file):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(id=' ',file=in_file)
    with open(out_file,'w') as handle:

        if mode == 'binder':
            for ind,frag_c in enumerate(cluster_coord_list([atom.get_coord() for atom in structure.get_atoms()])):
                # print "ATOM    100  LFC LFC X   0      -0.868 -10.433  16.680" \\this is the reference line
                handle.write("ATOM    100  LFC LFC X {:>3}    {:>8.3f}{:>8.3f}{:>8.3f}\n".format(ind,frag_c[0],frag_c[1],frag_c[2]))
        elif mode == 'peptide':
            residues = list(structure.get_residues())
            for residue in residues:
                bb_atom_c = gmc_coord_list([atom.get_coord() for atom in residue.get_atom() if atom.backbone])
                handle.write("ATOM    100  LBC {} {} {:>3}    {:>8.3f}{:>8.3f}{:>8.3f}\n".format(residue.get_resname(),
                                                                                        residue.get_parent().id,
                                                                                        residue.get_id()[1], bb_atom_c[0],
                                                                                        bb_atom_c[1], bb_atom_c[2]))
                if len([atom.get_coord() for atom in residue.get_atom() if not atom.backbone]) != 0:
                    sc_atom_c = gmc_coord_list([atom.get_coord() for atom in residue.get_atom() if not atom.backbone])
                    # print "ATOM    100  LFC LFC X   0      -0.868 -10.433  16.680"
                    handle.write("ATOM    100  LSC {} {} {:>3}    {:>8.3f}{:>8.3f}{:>8.3f}\n".format(residue.get_resname(),residue.get_parent().id,residue.get_id()[1], sc_atom_c[0],
                                                                                        sc_atom_c[1], sc_atom_c[2]))









    return


if __name__ == '__main__':
    import sys
    args = sys.argv[1:]

    help_info = "How to use this program:\n" \
                "python get_lfc.py [opt] [files]\n" \
                "opt: -p for peptide where every residue is clustered on it own\n" \
                "opt: -l for binder where all atoms are clustered together"
    try:
        if '-l' == args[0]:
            for files in args[1:-1]:
                get_lfc(files,'binder',args[-1])
        elif '-p' == args[0]:
            for files in args[1:-1]:
                get_lfc(files,'peptide',args[-1])
        else:
            print(help_info)
    except:
        print(help_info)




