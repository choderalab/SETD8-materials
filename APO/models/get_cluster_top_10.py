import numpy as np
import mdtraj as md

cluster_top_10 = np.load('cluster_top_10.npy')

for i in range(100):
    top_10 = cluster_top_10[i]
    frame_count = 0
    for frame in top_10:
        if frame_count == 0:
            traj = md.load('/cbio/jclab/home/rafal.wiewiora/repos/MSM_play/set8_apo_11707_11709_FINAL/data_cut_start_noH/' + frame[0])[int(frame[1])]
            frame_count += 1
        else:
            traj_ = md.load('/cbio/jclab/home/rafal.wiewiora/repos/MSM_play/set8_apo_11707_11709_FINAL/data_cut_start_noH/' + frame[0])[int(frame[1])]
            traj = traj.join(traj_)
    traj.save('cluster_centers_top10/%d.h5' % i)
    traj.save('cluster_centers_top10/%d.dcd' % i)
