import numpy as np
import glob
import mdtraj as md
from multiprocessing import Pool

x = glob.glob('/cbio/jclab/home/rafal.wiewiora/repos/MSM_play/set8_apo_11707_11709_FINAL/data_cut_start_noH_stride10/11709/*.h5')
num_threads = 16

def featurize_distances(fname):
    traj = md.load(fname)
    pairwise_distances, residue_pairs = md.compute_contacts(traj, scheme='closest-heavy')
    np.save('/cbio/jclab/home/rafal.wiewiora/repos/MSM_play/set8_apo_11707_11709_FINAL/data_cut_start_noH_stride10_featurized/dist/11709/%s' % fname[108:-3], pairwise_distances)
    np.save('/cbio/jclab/home/rafal.wiewiora/repos/MSM_play/set8_apo_11707_11709_FINAL/data_cut_start_noH_stride10_featurized/dist_res_pairs/11709/%s' % fname[108:-3], residue_pairs)
    
pool = Pool(num_threads)
pool.map(featurize_distances, x)
