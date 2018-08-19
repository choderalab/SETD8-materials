import numpy as np
import pyemma
pyemma.config.use_trajectory_lengths_cache = False
import glob
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('nclusters', type=int)
parser.add_argument('split', type=int)
args = parser.parse_args()
nclusters = args.nclusters
split = args.split

splits_train = np.load('/home/rafal.wiewiora/repos/MSM/set8_apo_11707_11709_FINAL/feat_choice_scoring_new_scheme_2/splits_train_dih.npy')
splits_test = np.load('/home/rafal.wiewiora/repos/MSM/set8_apo_11707_11709_FINAL/feat_choice_scoring_new_scheme_2/splits_test_dih.npy')

data = pyemma.coordinates.load(['tica_trajs/%d/%d.npy' % (split,x) for x in splits_train[split]])

clustering = pyemma.coordinates.cluster_kmeans(data=data, k=nclusters, max_iter=1000, stride=10)

dtrajs = clustering.dtrajs
np.save('dtrajs/%d_%d_train.npy' % (nclusters, split), dtrajs)

del data
del dtrajs
data = [np.load('tica_trajs/%d/%d.npy' % (split,x)) for x in splits_test[split]]
dtrajs = clustering.transform(data)
np.save('dtrajs/%d_%d_test.npy' % (nclusters, split), dtrajs)
