"""
Calculate the SASA

Cluster based on SASA

Use representative states for AlphaSpace


"""
import os

import mdtraj as md
import numpy as np
from msmbuilder.cluster import MiniBatchKMeans

from msmbuilder.featurizer import Featurizer
from msmbuilder.io import load_meta, preload_tops, save_trajs, save_generic, gather_metadata, save_meta, \
    NumberedRunsParser, load_trajs, load_generic, preload_top, backup

from msmbuilder.decomposition import PCA
from multiprocessing import Pool

# define Parser for metadata
from msmbuilder.io.sampling import sample_states

parser = NumberedRunsParser(
    traj_fmt="{run}.h5",
    step_ps=50,
)
meta = gather_metadata("trajs/*.h5", parser)
save_meta(meta)

#define custom SASA featurizer
class SasaFeaturizer(Featurizer):
    """Featurizer based on dihedral angles.

    This featurizer transforms a dataset containing MD trajectories into
    a vector dataset by representing each frame in each of the MD trajectories
    by a vector containing one or more of the backbone or side-chain dihedral
    angles, or the sin and cosine of these angles.

    """

    def __init__(self, mode='atom'):
        self.mode = mode

        known = {'residue', 'atom'}
        if not mode in known:
            raise ValueError('angles must be a subset of %s. you supplied %s' %
                             (str(known), str(mode)))

    def partial_transform(self, traj, indice=None):
        """Featurize an MD trajectory into a vector space via calculation
        of dihedral (torsion) angles

        Parameters
        ----------
        traj : mdtraj.Trajectory
            A molecular dynamics trajectory to featurize.

        Returns
        -------
        features : np.ndarray, dtype=float, shape=(n_samples, n_features)
            A featurized trajectory is a 2D array of shape
            `(length_of_trajectory x n_features)` where each `features[i]`
            vector is computed by applying the featurization function
            to the `i`th snapshot of the input trajectory.

        See Also
        --------
        transform : simultaneously featurize a collection of MD trajectories
        """

        sasa = md.shrake_rupley(traj, mode=self.mode)

        if indice is not None:
            if max(indice) > sasa.shape[0]:
                raise Exception()
            else:
                return sasa[:, np.array(indice)]
        else:
            return sasa


## Load
meta = load_meta()
tops = preload_tops(meta)
sasa_feat = SasaFeaturizer(mode='residue')


## Featurize logic
def feat(irow):
    i, row = irow
    traj = md.load(row['traj_fn'], top=tops[row['top_fn']])
    feat_traj = sasa_feat.partial_transform(traj)
    return i, feat_traj


## Do it in parallel
with Pool() as pool:
    ftrajs = dict(pool.imap_unordered(feat, meta.iterrows()))

## Save
save_trajs(ftrajs, 'ftrajs', meta)
save_generic(sasa_feat, 'featurizer.pickl')



## Load
pca = PCA()
# meta, ftrajs = load_trajs("ftrajs")

## Fit
pca.fit(ftrajs.values())

## Transform
ttrajs = {}
for k, v in ftrajs.items():
    ttrajs[k] = pca.partial_transform(v)

## Save
save_trajs(ttrajs, 'ttrajs', meta)
save_generic(pca, 'pca.pickl')


kmeans = MiniBatchKMeans(n_clusters=200)
kmeans.fit([traj for traj in ttrajs.values()])

## Transform
ktrajs = {}
for k, v in ttrajs.items():
    ktrajs[k] = kmeans.partial_transform(v)

## Save
print(kmeans.summarize())
save_trajs(ktrajs, 'ktrajs', meta)
save_generic(kmeans, 'kmeans.pickl')

## Sample

inds = sample_states(ttrajs,
                     kmeans.cluster_centers_,
                     k=1)

save_generic(inds, "cluster-sample-inds.pickl")

## Make trajectories
top = preload_top(meta)
out_folder = "cluster_samples"
backup(out_folder)
os.mkdir(out_folder)

traj = []
for state_i, state_inds in enumerate(inds):
    traj.append(md.join(
        md.load_frame(meta.loc[traj_i]['traj_fn'], index=frame_i, top=top)
        for traj_i, frame_i in state_inds
    ))

traj = md.join(traj)
traj.save("{}.pdb".format(out_folder))

