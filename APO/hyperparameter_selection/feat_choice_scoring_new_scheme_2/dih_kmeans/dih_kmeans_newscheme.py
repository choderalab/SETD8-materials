import numpy as np
import pyemma
pyemma.config.use_trajectory_lengths_cache = False
import glob
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('tica_type', type=str)
parser.add_argument('lag_time', type=int)
parser.add_argument('nclusters', type=int)
parser.add_argument('split', type=int)
args = parser.parse_args()
tica_type = args.tica_type
lag_time = args.lag_time
nclusters = args.nclusters
split = args.split

splits_train = np.load('/home/rafal.wiewiora/repos/MSM/set8_apo_11707_11709_FINAL/feat_choice_scoring_new_scheme_2/splits_train_dih.npy')
splits_test = np.load('/home/rafal.wiewiora/repos/MSM/set8_apo_11707_11709_FINAL/feat_choice_scoring_new_scheme_2/splits_test_dih.npy')

if tica_type == 'kinetic':

    data = pyemma.coordinates.load(['tica_trajs_kinetic_lag%d/%d/%d.npy' % (lag_time,split,x) for x in splits_train[split]])

    clustering = pyemma.coordinates.cluster_kmeans(data=data, k=nclusters, max_iter=1000)

    dtrajs = clustering.dtrajs
    np.save('tica_trajs_kinetic_lag%d/dtrajs/%d_%d_train.npy' % (lag_time, nclusters, split), dtrajs)

    del data
    del dtrajs
    data = [np.load('tica_trajs_kinetic_lag%d/%d/%d.npy' % (lag_time,split,x)) for x in splits_test[split]]
    dtrajs = clustering.transform(data)
    np.save('tica_trajs_kinetic_lag%d/dtrajs/%d_%d_test.npy' % (lag_time, nclusters, split), dtrajs)

elif tica_type == 'commute':

    data = pyemma.coordinates.load(['tica_trajs_commute_lag%d/%d/%d.npy' % (lag_time,split,x) for x in splits_train[split]])

    clustering = pyemma.coordinates.cluster_kmeans(data=data, k=nclusters, max_iter=1000)
    
    dtrajs = clustering.dtrajs
    np.save('tica_trajs_commute_lag%d/dtrajs/%d_%d_train.npy' % (lag_time, nclusters, split), dtrajs)

    del data
    del dtrajs
    data = [np.load('tica_trajs_commute_lag%d/%d/%d.npy' % (lag_time,split,x)) for x in splits_test[split]]
    dtrajs = clustering.transform(data)
    np.save('tica_trajs_commute_lag%d/dtrajs/%d_%d_test.npy' % (lag_time, nclusters, split), dtrajs)

else:
    raise Exception('Wrong TICA mapping type')
