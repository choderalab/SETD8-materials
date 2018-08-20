import numpy as np
import glob

trajs = glob.glob('/cbio/jclab/home/rafal.wiewiora/repos/MSM_play/set8_ligands_11708_11710/data_cut_start_proteinonly_noH_featurized/dih/*.npy')

for traj in trajs:
    name = traj.split('/')[-1]
    strided_traj = np.load(traj)[::10]
    np.save('/cbio/jclab/home/rafal.wiewiora/repos/MSM_play/set8_ligands_11708_11710/data_cut_start_proteinonly_noH_stride10_featurized/dih/%s' % name, strided_traj)
    
