import numpy as np
import glob
import mdtraj as md
from multiprocessing import Pool

x = glob.glob('/cbio/jclab/home/rafal.wiewiora/repos/MSM_play/set8_apo_11707_11709_FINAL/data_cut_start_noH_stride10_featurized/dist/11709/*.npy')
num_threads = 16

# following function copied from MSMBuilder
steepness = 5
center = 0.5
def _transform(distances):
        result = 1.0/(1+np.exp(steepness*(distances-center)))
        return result

def logistic_transform_distances(fname):
    traj = np.load(fname)
    traj_logistic = _transform(traj)
    np.save('/cbio/jclab/home/rafal.wiewiora/repos/MSM_play/set8_apo_11707_11709_FINAL/data_cut_start_noH_stride10_featurized/dist_soft/11709/%s' 
            % fname[124:], traj_logistic)

pool = Pool(num_threads)
pool.map(logistic_transform_distances, x)    
