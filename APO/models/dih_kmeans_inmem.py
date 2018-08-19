import numpy as np
import pyemma
pyemma.config.use_trajectory_lengths_cache = False
import glob
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('nclusters', type=int)
args = parser.parse_args()
nclusters = args.nclusters

#source = pyemma.coordinates.source(glob.glob('trajs_95/*.npy'))
data = pyemma.coordinates.load(glob.glob('trajs_95/*.npy'))

clustering = pyemma.coordinates.cluster_kmeans(data=data, k=nclusters, max_iter=1000, stride=10)

#stages = [source, clustering]
#pipeline = pyemma.coordinates.pipeline(stages, chunksize = 5000)

dtrajs = clustering.dtrajs
np.save('dtrajs_strided_fit_/%d.npy' % nclusters, dtrajs)

centers = clustering.clustercenters
np.save('cluster_centers_strided_fit_/%d.npy' % nclusters, centers)
