from .Snapshot import *
from .functions import _binCluster, _group
from .Configs import _COLOR_DICT, _COLOR_IDX
from .View import draw_sphere
from .Cluster import _DPocket


class Trajectory(list):
    beta_cluster_dist = Snapshot.beta_cluster_dist
    pocket_cluster_dist = Snapshot.pocket_cluster_dist

    def __init__(self, snapshots):
        super().__init__(snapshots)
        self.dpocket_beta_labels = None
        self.beta_mapping = None

    @property
    def beta_xyz(self):
        """
        First Generate a mapping of beta atom to snapshot index and per snapshot beta index,
        then

        Returns
        -------

        """

        self.map_betas()

        beta_xyz = []
        for snapshot in self:
            beta_xyz.extend(snapshot._beta_xyz)

        return np.array(beta_xyz)

    def map_betas(self):
        beta_mapping = []
        for i, snapshot in enumerate(self):
            beta_mapping.extend([[i, beta_idx] for beta_idx in range(snapshot.numBetas)])
        self.beta_mapping = np.array(beta_mapping)

    def map_beta(self, beta_index):
        """

        Parameters
        ----------
        beta_index: array or int
            total index of beta atoms

        Returns
        -------
        [snapshot_idx, snapshot_beta_idx]


        """
        if self.beta_mapping is None:
            self.map_betas()
        return self.beta_mapping[beta_index]

    @property
    def dpockets(self):
        if self.dpocket_beta_labels is None:
            self.gen_dpockets(clust_distance=4.7, bin_size=[1.8, 1.8, 1.8])
        for beta_indices in self.dpocket_beta_labels:
            yield _DPocket(self, beta_indices=beta_indices)

    def gen_dpockets(self, clust_distance, bin_size=None):
        beta_coords = self.beta_xyz
        if bin_size is None:
            beta_dpocket_labels = fcluster(linkage(self.beta_xyz, method='average'), clust_distance / 10,
                                           criterion='distance') - 1
        else:
            beta_dpocket_labels = _binCluster(coords=beta_coords, distance=clust_distance / 10, bin_size=bin_size / 10)
        self.dpocket_beta_labels = _group(beta_dpocket_labels)

    def screen_by_contact(self, ref_traj, cutoff=3.6):

        if hasattr(ref_traj, 'xyz'):
            ref_traj = ref_traj.xyz
        elif hasattr(ref_traj, 'shape'):
            pass
        else:
            raise TypeError('ref_traj has to be mdtraj or numpy array')

        if len(ref_traj.shape) == 2:
            ref_traj = np.expand_dims(ref_traj, axis=0)

        if ref_traj.shape[0] != len(self):
            print("provided trajectory has less snapshots({}) than in this universe({})".format(ref_traj.shape[0],
                                                                                                len(self)))
            for i in range(len(self)):
                self[i].calculateContact(ref_traj[0], cutoff=cutoff)
        else:

            for i in range(len(self)):
                self[i].calculateContact(coords=ref_traj[i], cutoff=cutoff)

    def draw_pocket(self, view, snapshot_idx: int = 0, contact_only=True, radius=1.0):
        """
        Draw and view the current snapshot in jupyter notebook.

        Parameters
        ----------
        view
        snapshot_idx

        Returns
        -------
        view
            nglview widget object
        """

        i = 0
        for pocket in self[snapshot_idx].pockets:
            if (contact_only and pocket.isContact) or not contact_only:
                draw_sphere(view, list(pocket.centroid * 10), radius=radius, color=_COLOR_DICT[_COLOR_IDX[i]])
                print("_Pocket with index of {} is {}".format(pocket.index, _COLOR_IDX[i]))
                i = (i + 1) % len(_COLOR_DICT)

    def draw_alpha(self, view, snapshot_idx: int = 0, contact_only=True, radius=0.2):
        i = 0
        for pocket in self[snapshot_idx].pockets:
            if (contact_only and pocket.isContact) or not contact_only:
                for alpha in pocket.alphas:
                    draw_sphere(view, list(alpha.centroid * 10), radius, color=_COLOR_DICT[_COLOR_IDX[i]])

                print("_Pocket with index of {} is {}".format(pocket.index, _COLOR_IDX[i]))
                i = (i + 1) % len(_COLOR_DICT)

    def draw_beta(self, view, snapshot_idx: int = 0, contact_only=True, radius=0.2):
        i = 0
        for pocket in self[snapshot_idx].pockets:
            if (contact_only and pocket.isContact) or not contact_only:
                for beta in pocket.betas:
                    draw_sphere(view, list(beta.centroid * 10), radius, color=_COLOR_DICT[_COLOR_IDX[i]])
                print("_Pocket with index of {} is {}".format(pocket.index, _COLOR_IDX[i]))
                i = (i + 1) % len(_COLOR_DICT)

    def draw_dpocket(self, view, dpockets=None, contact_only=True, radius=1.0, label='index'):
        """

        Parameters
        ----------
        view
        dpockets
        contact_only
        radius
        label: str
            types of label u wishes to display

        Returns
        -------

        """
        dpockets = self.dpockets if dpockets is None else dpockets

        i = 0
        for dpocket in dpockets:
            if (contact_only and dpocket.is_contact) or not contact_only:
                draw_sphere(view, list(dpocket.centroid * 10), radius=radius,
                            color=_COLOR_DICT[_COLOR_IDX[i % len(_COLOR_DICT)]])

                if label == 'index':
                    view.shape.add_text(list(dpocket.centroid * 10), _COLOR_DICT[_COLOR_IDX[i % len(_COLOR_DICT)]],
                                        4, "   {}".format(i))
                elif label == 'score':
                    view.shape.add_text(list(dpocket.centroid * 10), _COLOR_DICT[_COLOR_IDX[i % len(_COLOR_DICT)]],
                                        4, "   {}".format(str(np.round(np.average(dpocket.scores), 2))))
                elif label == 'space':
                    view.shape.add_text(list(dpocket.centroid * 10), _COLOR_DICT[_COLOR_IDX[i % len(_COLOR_DICT)]],
                                        4, "   {}".format(int(np.average(dpocket.spaces))))
                else:
                    print('Label {} not recognized as'.format(label), 'index', 'score', 'space')
                    view.shape.add_text(list(dpocket.centroid * 10), _COLOR_DICT[_COLOR_IDX[i % len(_COLOR_DICT)]],
                                        4, "   {}".format(i))
                i += 1
