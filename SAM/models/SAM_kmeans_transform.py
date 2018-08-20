import numpy as np
import glob
import scipy

cluster_centers = np.load('/cbio/jclab/home/rafal.wiewiora/repos/MSM_play/set8_apo_11707_11709_FINAL/dih_full_dset/cluster_centers_strided_fit/100.npy')

trajs = glob.glob('dih/SAM/trajs_95/*.npy')

for traj in trajs:
    dtraj, distances = scipy.cluster.vq.vq(np.load(traj), cluster_centers)
    np.save('dih/SAM/dtrajs/%s' % traj.split('/')[-1], dtraj)
    np.save('dih/SAM/distances/%s' % traj.split('/')[-1], distances)
