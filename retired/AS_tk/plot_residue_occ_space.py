"""
This program plot the dynamic residue occupancy of each residue. For this instance, 845-858 residues are used to locate
the involved pockets, while only 847-858 are used for AlphaSpace counting
The clustering result will be printed out as the program run

Output file info:

residue name - index

residue occ space average, occ space SD, percent time in contact

"""

import matplotlib.pyplot as plt

from PDB.alphaspace import AlphaSpace
from PDB.structure_clustering import *

"""
Options start
"""

# If you wish to use clustering
use_clustering = True
# how many clusters to show on the plot
max_cluster_index = 4
# how many clusters should there be in total
total_cluster_number = 5
#  which clustering method should you use


method = 'lig_bb_rmsd'
# occupied or unoccupied space to display in the plot
occ_unocc = 'occ'
output_log_filename = 'log.dat'

pocket_detection_residue_indices = list(range(88, 101))
pocket_plot_residue_indices = list(range(90, 101))
residue_to_visualize = list(range(95, 101))
color = ['red', 'grey', 'blue', 'pink', 'yellow', 'purple', 'brown']

width = 0.35  # the width of the bars
"""
options end
"""

if __name__ == '__main__':
    try:
        import sys, os

        target_system = sys.argv[1]
    except:
        target_system = os.getcwd()

folders = [os.path.join(target_system, d) for d in os.listdir(target_system) if
           os.path.isdir(os.path.join(target_system, d))]
directories = []
for folder in folders:
    directories += ([os.path.join(folder, f) for f in os.listdir(folder) if
                     os.path.isdir(os.path.join(folder, f))])

if use_clustering:
    distance_matrix = calculate_distance_matrix(folders=folders, method=method, )
    clusters, cluster_rep_ind = structure_cluster(folders, distance_matrix=distance_matrix,
                                                  cluster_number=total_cluster_number,
                                                  output_average_structure=True)
    cluster_rep_ind = [y - 1 for (x, y) in
                       sorted(zip(clusters, cluster_rep_ind), key=lambda item: len(item[0]), reverse=True)]
    clusters.sort(key=lambda item: len(item), reverse=True)
else:
    count = 0
    clusters = [[] for _ in folders]
    for i, folder in enumerate(folders):
        for subf in [os.path.join(folder, f) for f in os.listdir(folder) if os.path.isdir(os.path.join(folder, f))]:
            clusters[i].append(count)
            count += 1
    cluster_rep_ind = list(range(len(clusters)))
    for i in range(len(clusters)):
        cluster_rep_ind[i] = structure_cluster(folders=[folders[i]],
                                               distance_matrix=calculate_distance_matrix(folders=[folders[i], ],
                                                                                         method=method),
                                               cluster_number=1, output_average_structure=True)[1][0]

res_ids = set()
for folder in directories:
    with open(os.path.join(folder, 'table_cntctscore_new.dat')) as handle:
        for line in handle.readlines():
            res_id = line.split()[0]
            if res_id not in ['resID', 'FACE'] and int(res_id[:3]) in pocket_plot_residue_indices:
                res_ids.add(res_id)

res_ids = sorted([id for id in res_ids])
ind = np.arange(len(res_ids))

fig, ax = plt.subplots()
ax2 = ax.twinx()

cluster_bar = []
cluster_scatter = []

with open(output_log_filename, 'w') as handle:
    if use_clustering:
        handle.write("cluster index:\n")
    else:
        handle.write("folder index\n")
    for c_i, cluster in enumerate(clusters):
        handle.write(str(cluster_rep_ind[c_i]) + " : ")
        for i in cluster:
            handle.write(str(i) + "   ")
        handle.write("\n")
    handle.write("\n")

y_lim = 0.0
for cluster_index in range(0, min(max_cluster_index, len(clusters))):
    with open(output_log_filename, 'a') as handle:
        handle.write("cluster " + str(cluster_index + 1) + ':\n')
    residue_occurance_count = {}
    residue_unoccupied_count = {}
    residue_space = {}
    for residue in res_ids:
        residue_space[residue] = []
        residue_occurance_count[residue] = 0
        residue_unoccupied_count[residue] = 0
    for i, dir in enumerate([directories[d] for d in clusters[cluster_index]]):
        current_alphaspace = AlphaSpace(folder_dir=dir)
        ligand_atoms = current_alphaspace.get_atoms(structure='binder')
        ligand = current_alphaspace.get_structure(structure='binder')

        # create all pocket list
        pocket_point_list = []
        residue_linked_pocket = {}
        for pocket in current_alphaspace.contact_pocket:
            for pocket_point in pocket.get_points(type='AAC'):
                pocket_point_list.append(pocket_point)

        # calculate the distance matrix
        ppdist_mat = cdist([pocket_point.coord for pocket_point in pocket_point_list],
                           [atom.get_coord() for atom in ligand_atoms])

        # select the association ones: close to 845-858
        selected_pocket_point_list = []
        for i, pocket_point_dist in enumerate(ppdist_mat):
            linked_residue_index = ligand_atoms[pocket_point_dist.argmin()].get_parent().get_id()[1]
            if linked_residue_index in pocket_detection_residue_indices:
                selected_pocket_point_list.append(pocket_point_list[i])
        # recalculate the distance matrix for the selected pockets
        sel_ppdist_mat = cdist([pocket_point.coord for pocket_point in selected_pocket_point_list],
                               [atom.get_coord() for atom in ligand_atoms if
                                atom.get_parent().get_id()[1] in pocket_plot_residue_indices])

        for i, pocket_point_dist in enumerate(sel_ppdist_mat):
            linked_residue_index = \
                [atom for atom in ligand_atoms if atom.get_parent().get_id()[1] in pocket_plot_residue_indices][
                    pocket_point_dist.argmin()].get_parent().get_id()[1]
            if linked_residue_index in residue_linked_pocket:
                residue_linked_pocket[linked_residue_index].append(selected_pocket_point_list[i])
            else:
                residue_linked_pocket[linked_residue_index] = [selected_pocket_point_list[i]]

        for i in residue_space:

            current_residue = [residue for residue in ligand.get_residues() if residue.get_id()[1] == int(i[:3])][0]
            occupancy = [0.0, 0.0]
            temp = 0.0
            if int(i[:3]) in residue_linked_pocket:
                for pocket_point in residue_linked_pocket[int(i[:3])]:
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
    res_ids = [res_id for res_id in res_ids if int(res_id[:3]) in residue_to_visualize]

    residue_average_occ_space = [np.average(residue_space[residue], axis=0)[0] for residue in res_ids]
    residue_average_unocc_space = [np.average(residue_space[residue], axis=0)[1] for residue in res_ids]
    residue_occ_sd = [np.std(residue_space[residue], axis=0)[0] for residue in res_ids]
    residue_unocc_sd = [np.std(residue_space[residue], axis=0)[1] for residue in res_ids]
    residue_contact_time = [residue_occurance_count[residue] / len(clusters[cluster_index])
                            for residue in res_ids]
    ind = np.arange(len(res_ids))
    if occ_unocc == 'occ':
        cluster_bar.append(ax.bar(ind + cluster_index * width * 2 / max_cluster_index,
                              residue_average_occ_space,
                              width * 2 / max_cluster_index,
                              color=color[cluster_index], ecolor='black', edgecolor=color[cluster_index],
                              yerr=residue_occ_sd))


        if y_lim < max(residue_average_occ_space):
            y_lim = max(residue_average_occ_space)
        cluster_scatter.append(
            ax2.scatter(ind + cluster_index * width * 2 / max_cluster_index + width / max_cluster_index, residue_contact_time,
                        s=100, c=color[cluster_index], alpha=1.0))

    else:
        cluster_bar.append(ax.bar(ind + cluster_index * width * 2 / max_cluster_index,
                                  residue_average_unocc_space,
                                  width * 2 / max_cluster_index,
                                  color=color[cluster_index], ecolor='black', edgecolor=color[cluster_index],
                                  yerr=residue_unocc_sd))
        if y_lim < max(residue_average_unocc_space):
            y_lim = max(residue_average_unocc_space)
        cluster_scatter.append(ax2.scatter(ind + cluster_index * width * 2 / max_cluster_index + width / 4,
                                           [residue_unoccupied_count[residue] / len(clusters[cluster_index])
                                            for residue in res_ids],
                                           s=100, c=color[cluster_index], alpha=1.0))

    with open(output_log_filename, 'a') as handle:
        handle.write("residue average_occ_space occ_sd average_unocc_space unocc_sd percent_time_in_contact\n")
        for i in range(len(residue_average_occ_space)):
            handle.write("{} {:3d} {:3d} {:3d} {:3d} {:3d}%\n".format(res_ids[i], int(residue_average_occ_space[i]),
                                                                      int(residue_occ_sd[i]),
                                                                      int(residue_average_unocc_space[i]),
                                                                      int(residue_unocc_sd[i]),
                                                                      int(residue_contact_time[i] * 100)))
        handle.write('\n')

ax.set_xticks(ind + width)
ax.set_ylabel('average_space_' + occ_unocc, fontsize=22)
ax.set_ylim((0, y_lim + 10))
plt.xlim([0, ind.size])
ax2.set_ylim((0, 1.1))
plt.legend((cluster_bar), ['Cluster ' + str(i + 1) + ' ' + str(len(clusters[i])) for i in
                           range(min(max_cluster_index, len(clusters)))])
# plt.legend((cluster_bar),['Cluster '+ str(i+1)+' '+ str(len(clusters[i])) for i in range(1,min(max_cluster_index,len(clusters)))])
if occ_unocc == 'occ':
    ax2.set_ylabel('Percent time in contact', fontsize=22)
else:
    ax2.set_ylabel('Percent time in association with unoccupied space', fontsize=22)
ax.set_xticklabels(res_ids)
plt.title(method + '   ' + folders[0][:-4])
plt.show()
