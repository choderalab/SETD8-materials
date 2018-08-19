import pyemma
import glob
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('split', type=int)
args = parser.parse_args()
split = args.split

splits_train = np.load('/home/rafal.wiewiora/repos/MSM/set8_apo_11707_11709_FINAL/feat_choice_scoring_new_scheme_2/splits_train_dist.npy')
#splits_test = np.load('/cbio/jclab/home/rafal.wiewiora/repos/MSM_play/set8_apo_11707_11709_FINAL/feat_choice_scoring_new_scheme/splits_test_dist.npy')

source = pyemma.coordinates.source(list(np.array(glob.glob('/data/chodera/rafal.wiewiora/set8_apo_11707_11709/data_cut_start_noH_stride10_featurized/dist_cross/*/*.npy'))[list(splits_train[split])]))

tica_commute = pyemma.coordinates.tica(lag=1, kinetic_map=False, commute_map=True, var_cutoff=1)
stages = [source, tica_commute]
pipeline = pyemma.coordinates.pipeline(stages, chunksize = 1000)
dim_commute = np.argwhere(np.cumsum(tica_commute.timescales)/np.sum(tica_commute.timescales) > 0.95)[0,0] + 1

for i, traj in enumerate(glob.glob('/data/chodera/rafal.wiewiora/set8_apo_11707_11709/data_cut_start_noH_featurized/data_cut_start_noH_featurized/dist/*/*.npy')):
    traj = np.load(traj)
    traj = tica_commute.transform(traj)
    traj = traj[:,:dim_commute]
    np.save('tica_trajs/%d/%d.npy' % (split,i), traj)
