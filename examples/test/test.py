from alphaspace import *
import sys




if __name__ == '__main__':

    universe = get_path(sys.argv[1])

    beta_xyz = universe.beta_xyz()

    bin_size = [1.8, 1.8, 1.8]

    beta_pocket_labels = bin_cluster(beta_xyz * 10, bin_size=bin_size)

    # remap the label to each snapshot

    for beta_idx, pocket_idx in enumerate(beta_pocket_labels):
        snapshot_idx, snapshot_beta_idx = universe.map_beta(beta_idx)
        probe_scores = universe[snapshot_idx].beta_scores[snapshot_beta_idx]

        best_probe_element = probe_element_list[np.argmin(probe_scores)]
        best_probe_score = probe_scores[np.argmin(probe_scores)]

        print("Beta Atom {} in snapshot {} belongs to pocket {} with {} of score {}".format(snapshot_beta_idx,
                                                                                            snapshot_idx, pocket_idx,
                                                                                            best_probe_element,
                                                                                            best_probe_score))
