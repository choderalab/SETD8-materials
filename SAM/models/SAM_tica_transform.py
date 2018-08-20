# SAM is: project 11710 RUN0 (TS_SAM) and RUN2 (SAM_SAM)

import pyemma
import glob
import numpy as np
import mdtraj as md

lag_time = 1

# remove the unfolding trajectory - number 3959
source = pyemma.coordinates.source([x for x in glob.glob('/cbio/jclab/home/rafal.wiewiora/repos/MSM_play/set8_apo_11707_11709_FINAL/data_cut_start_noH_stride10_featurized/dih/*.npy') if '3959' not in x])

tica_commute = pyemma.coordinates.tica(lag=lag_time, kinetic_map=False, commute_map=True, var_cutoff=1)
stages = [source, tica_commute]
pipeline = pyemma.coordinates.pipeline(stages, chunksize = 1000)

trajs = glob.glob('data_cut_start_proteinonly_noH_featurized/dih/SAM/*.npy')

for traj in trajs:
    name = traj.split('/')[-1]
    traj_ = tica_commute.transform(np.load(traj))
    np.save('dih/SAM/trajs/%s' % name, traj_)
    traj_ = traj_[:, :467]
    np.save('dih/SAM/trajs_95/%s' % name, traj_)
    traj_ = traj_[:, :10]
    np.save('dih/SAM/trajs_10/%s' % name, traj_)
