import pickle
import os
import shutil
import sys

import matplotlib.pyplot as plt

from ..PDB.alphaspace import AlphaSpace,find_second_last
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
        if folder[:len(option.alphaspace_folder_prefix)] == option.alphaspace_folder_prefix and os.path.isdir(
                os.path.join(cwd, folder)):
            as_master_dirs.append(os.path.join(cwd, folder))
else:
    for folder in os.listdir(cwd):
        if os.path.isdir(os.path.join(cwd, folder)):
            as_master_dirs.append(os.path.join(cwd, folder))
for as_master_dir in as_master_dirs:
    for as_sub_dir in os.listdir(as_master_dir):
        if os.path.isdir(os.path.join(as_master_dir, as_sub_dir)):
            as_sub_dirs.append(os.path.join(as_master_dir, as_sub_dir))

if option.use_clustering:
    if "dist_mat_catche.dat" in os.listdir(cld):
        print("Existing distance matrix found")
        with open(os.path.join(cld, 'dist_mat_catche.dat'), 'r') as handle:
            dist_mat_cache_file = pickle.load(handle)
    else:
        print("No distance matrix found, creating new one")
        dist_mat_cache_file = calculate_distance_matrix(folders=as_master_dirs, method=option.clustering_method)
        with open(os.path.join(cld, 'dist_mat_catche.dat'), 'w') as handle:
            pickle.dump(dist_mat_cache_file, handle)
    clusters, cluster_rep_ind = structure_cluster(as_master_dirs, distance_matrix=dist_mat_cache_file,
                                                  cluster_number=option.total_cluster_number,
                                                  output_average_structure=True)
    cluster_rep_ind = [y - 1 for (x, y) in
                       sorted(zip(clusters, cluster_rep_ind), key=lambda item: len(item[0]), reverse=True)]
    clusters.sort(key=lambda item: len(item), reverse=True)
    clusters = clusters[:option.total_cluster_number]
    cluster_rep_ind = cluster_rep_ind[:option.total_cluster_number]

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
                                                           distance_matrix=calculate_distance_matrix(
                                                               folders=[as_master_dirs[i], ],
                                                               method=option.clustering_method),
                                                           cluster_number=1, output_average_structure=True)[1][0]]

res_ids = option.residue_to_visualize
ind = np.arange(len(res_ids))

fig, ax = plt.subplots()
residue_space_log = {}
occ_bar = []
unocc_bar = []
y_lim = option.y_lim
res_name_dict = {}


with open(os.path.join(cld,option.output_log_filename), 'w') as handle:
    if option.use_clustering:
        handle.write("cluster index:\n")
    else:
        handle.write("folder index\n")
    for c_i, cluster in enumerate(clusters):
        handle.write(str(cluster_rep_ind[c_i]) + " : ")
        for i in cluster:
            handle.write(str(i) + "   ")
        handle.write("\n")
    handle.write("\n")


for cluster_index in range(0, min(option.max_cluster_index, len(clusters))):
    with open(os.path.join(cld,option.output_log_filename), 'a') as handle:
        handle.write("cluster " + str(cluster_index + 1) + ':\n')
    print(len(clusters[cluster_index]))
    residue_occurance_count = {}
    residue_unoccupied_count = {}
    residue_space = {}
    for residue in res_ids:
        residue_space[residue] = []
        residue_occurance_count[residue] = 0
        residue_unoccupied_count[residue] = 0
    for directory in [as_sub_dirs[d] for d in clusters[cluster_index]]:
        current_alphaspace = AlphaSpace(folder_dir=directory)
        ligand_atoms = current_alphaspace.get_atoms(structure='binder')
        ligand = current_alphaspace.get_structure(structure='binder')
        for residue in ligand.get_residues():
            res_name_dict[residue.get_id()[1]] = residue.get_resname()

        # create all pocket list
        pocket_point_list = []
        residue_linked_pocket = {}
        for pocket in current_alphaspace.contact_pocket:
            for pocket_point in pocket.get_points(type='AAC'):
                pocket_point_list.append(pocket_point)

        # calculate the distance matrix
        ppdist_mat = cdist([pocket_point.coord for pocket_point in pocket_point_list],
                           [atom.get_coord() for atom in ligand_atoms])

        # select the association ones:
        selected_pocket_point_list = []
        for i, pocket_point_dist in enumerate(ppdist_mat):
            linked_residue_index = ligand_atoms[pocket_point_dist.argmin()].get_parent().get_id()[1]
            if linked_residue_index in option.pocket_detection_residue_indices:
                selected_pocket_point_list.append(pocket_point_list[i])
        # recalculate the distance matrix for the selected pockets
        sel_ppdist_mat = cdist([pocket_point.coord for pocket_point in selected_pocket_point_list],
                               [atom.get_coord() for atom in ligand_atoms if
                                atom.get_parent().get_id()[1] in option.pocket_plot_residue_indices])

        for i, pocket_point_dist in enumerate(sel_ppdist_mat):
            linked_residue_index = \
                [atom for atom in ligand_atoms if atom.get_parent().get_id()[1] in option.pocket_plot_residue_indices][
                    pocket_point_dist.argmin()].get_parent().get_id()[1]
            if linked_residue_index in residue_linked_pocket:
                residue_linked_pocket[linked_residue_index].append(selected_pocket_point_list[i])
            else:
                residue_linked_pocket[linked_residue_index] = [selected_pocket_point_list[i]]

        for i in residue_space:
            current_residue = [residue for residue in ligand.get_residues() if residue.get_id()[1] == i][0]
            occupancy = [0.0, 0.0]
            temp = 0.0
            if i in residue_linked_pocket:
                for pocket_point in residue_linked_pocket[i]:
                    if min(cdist([pocket_point.coord], [atom.get_coord() for atom in current_residue])[0]) <= 1.6:
                        occupancy[0] = occupancy[0] + pocket_point.polar_score + pocket_point.nonpolar_score
                    else:
                        occupancy[1] = occupancy[1] + pocket_point.polar_score + pocket_point.nonpolar_score
                        temp = 1.0
            if occupancy[0] > 0.0:
                residue_occurance_count[i] += 1.0
            residue_unoccupied_count[i] += temp
            residue_space[i].append(occupancy)

    # pick residue to visualize
    res_ids = [res_id for res_id in res_ids if res_id in option.residue_to_visualize]

    residue_average_occ_space = [np.average(residue_space[residue], axis=0)[0] for residue in res_ids]
    residue_average_unocc_space = [np.average(residue_space[residue], axis=0)[1] for residue in res_ids]
    residue_occ_sd = [np.std(residue_space[residue], axis=0)[0] for residue in res_ids]
    residue_unocc_sd = [np.std(residue_space[residue], axis=0)[1] for residue in res_ids]
    residue_contact_time = [residue_occurance_count[residue] / len(clusters[cluster_index])
                            for residue in res_ids]
    ind = np.arange(len(res_ids))


    # save resdiue space log to file
    residue_space_log_folder = os.path.join(cld, "res_space_log_"+str(cluster_index))
    if "res_space_log_"+str(cluster_index) in os.listdir(cld):
        shutil.rmtree(residue_space_log_folder)
    os.mkdir(residue_space_log_folder)
    for residue in res_ids:
        with open(os.path.join(residue_space_log_folder,str(residue)+".dat"),'w') as handle:
            handle.write("cluster_ind occ_space unocc_space directory\n")
            for dir_ind, dir in enumerate([as_sub_dirs[d] for d in clusters[cluster_index]]):
                handle.write("{}\t{:5}\t{:5}\t{}\n".format(cluster_index,residue_space[residue][dir_ind][0],residue_space[residue][dir_ind][1],dir[find_second_last(dir,'/'):]))
    with open(os.path.join(residue_space_log_folder,"FACE.dat"), 'w') as handle:
        handle.write("cluster_ind occ_space unocc_space directory\n")
        for dir_ind, dir in enumerate([as_sub_dirs[d] for d in clusters[cluster_index]]):
            occ_sum,unocc_sum = np.sum([residue_space[residue][dir_ind] for residue in res_ids],axis=0)
            handle.write("{}\t{:5}\t{:5}\t{}\n".format(cluster_index, occ_sum,
                                                       unocc_sum,
                                                       dir[find_second_last(dir, '/'):]))


    if option.error_bar:

        occ_bar.append(ax.bar(ind + cluster_index * option.width * 2 / option.max_cluster_index,
                              residue_average_occ_space,
                              option.width * 2 / option.max_cluster_index,
                              color=option.color[cluster_index], ecolor='black', edgecolor=option.color[cluster_index],
                              yerr=residue_occ_sd))
        unocc_bar.append(ax.bar(ind + cluster_index * option.width * 2 / option.max_cluster_index,
                                residue_average_unocc_space,
                                option.width * 2 / option.max_cluster_index,
                                color=option.color[cluster_index], ecolor='black', bottom=residue_average_occ_space,
                                edgecolor=option.color[cluster_index],
                                yerr=residue_unocc_sd, alpha=0.5))
    else:
        occ_bar.append(ax.bar(ind + cluster_index * option.width * 2 / option.max_cluster_index,
                              residue_average_occ_space,
                              option.width * 2 / option.max_cluster_index,
                              color=option.color[cluster_index], ecolor='black',
                              edgecolor=option.color[cluster_index], ))
        unocc_bar.append(ax.bar(ind + cluster_index * option.width * 2 / option.max_cluster_index,
                                residue_average_unocc_space,
                                option.width * 2 / option.max_cluster_index,
                                color=option.color[cluster_index], ecolor='black', bottom=residue_average_occ_space,
                                edgecolor=option.color[cluster_index], alpha=0.5))


    with open(os.path.join(cld,option.output_log_filename), 'a') as handle:
        handle.write("residue average_occ_space occ_sd average_unocc_space unocc_sd percent_time_in_contact\n")
        for i in range(len(residue_average_occ_space)):
            handle.write("{} {:3d} {:3d} {:3d} {:3d} {:3d}%\n".format(res_ids[i], int(residue_average_occ_space[i]),
                                                                      int(residue_occ_sd[i]),
                                                                      int(residue_average_unocc_space[i]),
                                                                      int(residue_unocc_sd[i]),
                                                                      int(residue_contact_time[i] * 100)))
        handle.write('\n')

ax.set_xticks(ind + option.width)
ax.set_xticklabels([str(res_id) + "-" + res_name_dict[res_id] for res_id in option.residue_to_visualize])

ax.set_xticklabels = [i * 20 for i in range(10)]
ax.set_yticks(np.array([i * 20 for i in range(10)]))

for t in ax.get_xticklabels():
    t.set_fontsize(18)
for t in ax.get_yticklabels():
    t.set_fontsize(18)

ax.set_ylabel("occupied/unoccupied alpha-space", fontsize=28)
ax.set_ylim((0, y_lim))
plt.xlim([0, ind.size + option.width])

# plt.legend((occ_bar[0], unocc_bar[0], occ_bar[1], unocc_bar[1]),
#        ["binding mode 1 (occupied)", "binding mode 2 (occupied)", "binding mode 1 (unoccupied)",
#         "binding mode 2 (unoccupied)"], prop={'size': 18})

fig.set_dpi(300)
ax.yaxis.grid()
fig.set_size_inches(22, 10)
plt.savefig(os.path.join(cld,"plot.png"), bbox_inches='tight', )
