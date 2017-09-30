""""""

import pickle

from ..PDB.alphaspace import *
from structure_cluster import *

cwd = os.getcwd()
sys.path.append(cwd)
try:
    option_file_module_name = sys.argv[1][:-3]
    option = __import__(option_file_module_name)
    print("Using custom option "+ option_file_module_name)
except:
    if "option.py" not in os.listdir(cwd):
        print("Option file not found, using default options.")
        write_option_file(os.path.join(cwd, 'option.py'), "option_backup.py")
    # noinspection PyUnresolvedReferences
    import option as option
    option_file_module_name = "option"


cld = os.path.join(cwd, option_file_module_name)
try:
    os.mkdir(cld)
except:
    pass





as_master_dirs = []
as_sub_dirs = []
if option.alphaspace_folder_prefix is not None and len(option.alphaspace_folder_prefix) != 0:
    for folder in os.listdir(cwd):
        if folder[:len(option.alphaspace_folder_prefix)] == option.alphaspace_folder_prefix and os.path.isdir(os.path.join(cwd,folder)):
            as_master_dirs.append(os.path.join(cwd,folder))
else:
    for folder in os.listdir(cwd):
        if os.path.isdir(os.path.join(cwd,folder)):
            as_master_dirs.append(os.path.join(cwd,folder))
for as_master_dir in as_master_dirs:
    for as_sub_dir in os.listdir(as_master_dir):
        if os.path.isdir(os.path.join(as_master_dir,as_sub_dir)):
            as_sub_dirs.append(os.path.join(as_master_dir,as_sub_dir))



dist_mat_cache_file = None
if option.use_clustering:
    if "dist_mat_catche.dat" in os.listdir(cld):
        print("Existing distance matrix found")
        with open(os.path.join(cld, 'dist_mat_catche.dat'), 'r') as handle:
            dist_mat_cache_file = pickle.load(handle)
    else:
        print("No distance matrix found, creating new one")
        dist_mat_cache_file = calculate_distance_matrix(folders=as_master_dirs, method=option.clustering_method)
        with open(os.path.join(cld,'dist_mat_catche.dat'),'w') as handle:
            pickle.dump(dist_mat_cache_file,handle)
    clusters, cluster_rep_ind = structure_cluster(as_master_dirs, distance_matrix=dist_mat_cache_file,
                                                  cluster_number=option.cluster_number,
                                                  output_average_structure=True)
    cluster_rep_ind = [y - 1 for (x, y) in
                       sorted(zip(clusters, cluster_rep_ind), key=lambda item: len(item[0]), reverse=True)]
    clusters.sort(key=lambda item: len(item), reverse=True)
    clusters = clusters[:option.cluster_number]
    cluster_rep_ind = cluster_rep_ind[:option.cluster_number]

else:
    count = 0
    clusters = [[] for _ in as_master_dirs]
    for i, folder in enumerate(as_master_dirs):
        for subf in [os.path.join(folder, f) for f in os.listdir(folder) if os.path.isdir(os.path.join(folder, f))]:
            clusters[i].append(count)
            count += 1
    cluster_rep_ind = [0 for _ in range(len(clusters))]
    for i in range(len(clusters)):
        cluster_rep_ind[i] = clusters[i][structure_cluster(folders=[as_master_dirs[i]],
                                               distance_matrix=calculate_distance_matrix(folders=[as_master_dirs[i], ],
                                                                                         method=option.clustering_method),
                                               cluster_number=1, output_average_structure=True)[1][0]]

for cluster_ind, cluster in enumerate(clusters):
    residue_linked_pocket = {}
    for folder in [as_sub_dirs[i] for i in cluster]:
        #  load alpha space into class object
        current_alphaspace = AlphaSpace(folder_dir=folder)
        ligand_atoms = current_alphaspace.get_atoms(structure='binder')
        ligand = current_alphaspace.get_structure(structure='binder')

        # put all ACC in a list
        pocket_point_list = []
        for pocket in current_alphaspace.contact_pocket:
            for pp in pocket.get_points(type='ACC'):
                pocket_point_list.append(pp)
                # print pocket_point.text_line
                # print pocket_point.text_line[:31]

        # link every ACC with nearest residue
        dist_matrix = cdist([p.get_coord() for p in pocket_point_list], [atom.get_coord() for atom in ligand_atoms])
        for j, i in enumerate(np.argmin(dist_matrix, axis=1)):
            pocket_point_list[j].update_occupy(linked_residue_index=ligand_atoms[i].get_parent().get_id()[1])
        for res_id in sorted([residue.get_id()[1] for residue in ligand.get_residues()]):
            if res_id in residue_linked_pocket:
                residue_linked_pocket[res_id] += [ACC for ACC in pocket_point_list if ACC.residue_index == res_id]
            else:
                residue_linked_pocket[res_id] = [ACC for ACC in pocket_point_list if ACC.residue_index == res_id]

    if len(option.residue_list) != 0:
        for res_id in [_ for _ in residue_linked_pocket]:
            if res_id not in option.residue_list:
                residue_linked_pocket.pop(res_id,None)

    with open(os.path.join(cld,"cluster_dcc_" + str(cluster_ind)+".pdb"), 'w') as handle:
        with open(os.path.join(cld,"cluster_log_" + str(cluster_ind)+".txt"), 'w') as log_handle:
            if option.use_clustering:
                for ind in cluster:
                    log_handle.write(str(ind) + "  ")
                log_handle.write('\n')
            log_handle.write(
                'representative structure is: ' + str(cluster_rep_ind[cluster_ind]) + "\n")
            log_handle.write(as_sub_dirs[cluster_rep_ind[cluster_ind]] + '\n')
            log_handle.write('residue\tspread\tpopulation\n')
            for res_id in residue_linked_pocket:
                if len(residue_linked_pocket[res_id]) != 0:
                    for ACC in residue_linked_pocket[res_id]:
                        handle.write(ACC.text_line)
                    x, y, z = np.average([ACC.get_coord() for ACC in residue_linked_pocket[res_id]], axis=0)

                    if option.spread_mode == "AVERAGE":
                        log_handle.write(space_fill(res_id) + "\t{:3.1f}".format(np.average(
                            [np.sum((ACC.get_coord() - np.array([x, y, z])) ** 2) ** 0.5 for ACC in
                             residue_linked_pocket[res_id]]))+ "\t{:3d}\n".format(len(residue_linked_pocket[res_id])))

                    elif option.spread_mode == 'RMSD':
                        log_handle.write(space_fill(res_id) + "\t{:3.1f}\n".format(np.sqrt(np.sum(np.square(
                            [(ACC.get_coord() - np.array([x, y, z])) for ACC in residue_linked_pocket[res_id]]))) / len(
                            residue_linked_pocket[res_id]))+ "\t{:3d}\n".format(len(residue_linked_pocket[res_id])))
                    DCC = pocket_point(point_type='DCC', coordinate=[x, y, z])
                    DCC.text_line = "ATOM      1  DCC DCC X   0    " + space_fill("%.3f" % x, 8) + space_fill(
                        "%.3f" % y, 8) + space_fill("%.3f" % z, 8)
                    DCC.update_occupy(linked_residue_index=res_id)
                    handle.write(DCC.text_line)
                    handle.write('\n')
