from sklearn.datasets import make_blobs
import pandas as pd
import hdbscan

blobs, labels = make_blobs(n_samples=2000, n_features=3, centers=10)

print(blobs.shape)

clusterer = hdbscan.HDBSCAN(metric='euclidean')

clusterer.fit(blobs)

print(clusterer.labels_.max())

from sklearn.metrics.pairwise import pairwise_distances

distance_matrix = pairwise_distances(blobs, metric='euclidean')
clusterer = hdbscan.HDBSCAN(metric='precomputed')
clusterer.fit(distance_matrix)

print(clusterer.labels_.max())
