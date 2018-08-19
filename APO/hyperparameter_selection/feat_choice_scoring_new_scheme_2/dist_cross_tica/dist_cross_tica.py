import pyemma
import glob
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('tica_type', type=str)
parser.add_argument('lag_time', type=int)
parser.add_argument('split', type=int)
args = parser.parse_args()
tica_type = args.tica_type
lag_time = args.lag_time
split = args.split

splits_train = np.load('/home/rafal.wiewiora/repos/MSM/set8_apo_11707_11709_FINAL/feat_choice_scoring_new_scheme_2/splits_train_dist.npy')
#splits_test = np.load('/cbio/jclab/home/rafal.wiewiora/repos/MSM_play/set8_apo_11707_11709_FINAL/feat_choice_scoring_new_scheme/splits_test_dist.npy')

source = pyemma.coordinates.source(list(np.array(glob.glob('/data/chodera/rafal.wiewiora/set8_apo_11707_11709/data_cut_start_noH_stride10_featurized/dist_cross/*/*.npy'))[list(splits_train[split])]))

# kinetic map
if tica_type == 'kinetic':

    tica_kinetic = pyemma.coordinates.tica(lag=lag_time, kinetic_map=True, var_cutoff=1)
    stages = [source, tica_kinetic]
    pipeline = pyemma.coordinates.pipeline(stages, chunksize = 1000)

    for i, traj in enumerate(glob.glob('/data/chodera/rafal.wiewiora/set8_apo_11707_11709/data_cut_start_noH_stride10_featurized/dist_cross/*/*.npy')):
        traj = np.load(traj)
        traj = tica_kinetic.transform(traj)
        np.save('tica_trajs_kinetic_lag%d/%d/%d.npy' % (lag_time,split,i), traj)  

elif tica_type == 'commute':

    tica_commute = pyemma.coordinates.tica(lag=lag_time, kinetic_map=False, commute_map=True, var_cutoff=1)
    stages = [source, tica_commute]
    pipeline = pyemma.coordinates.pipeline(stages, chunksize = 1000)

    for i, traj in enumerate(glob.glob('/data/chodera/rafal.wiewiora/set8_apo_11707_11709/data_cut_start_noH_stride10_featurized/dist_cross/*/*.npy')):
        traj = np.load(traj)
        traj = tica_commute.transform(traj)
        np.save('tica_trajs_commute_lag%d/%d/%d.npy' % (lag_time,split,i), traj)

else:
    raise Exception('Wrong TICA mapping type')        
