import glob
import sklearn
import sklearn.metrics
import numpy as np

trajs = glob.glob('trajs_95/*.npy')
cluster_centers = np.load('cluster_centers_strided_fit/100.npy')

for traj in trajs:
    name = traj.split('/')[-1]
    traj_ = np.load(traj)
    distances = sklearn.metrics.pairwise.pairwise_distances(traj_, Y=cluster_centers, n_jobs=16)
    np.save('cluster_distances/%s' % name, distances)

# now we're going to limit the distances to the minimum of the frame - meaning only the distance to its
# cluster center
for traj in trajs:
    name = traj.split('/')[-1]
    distances = np.load('cluster_distances/%s' % name)
    distances_min = np.amin(distances, axis=1)
    np.save('cluster_distances_min/%s' % name, distances_min)
    
