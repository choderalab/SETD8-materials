import pyemma
import mdtraj as md
import glob
import numpy as np
from multiprocessing import Pool

fnames = glob.glob('data_cut_start_noH_stride10_featurized/dist/*/*.npy')
res_pairs_fnames = glob.glob('data_cut_start_noH_stride10_featurized/dist_res_pairs/*/*.npy')

# first check that the residue pairs are the same order in all trajectories
x = np.load(res_pairs_fnames[0])
for i in range(len(res_pairs_fnames)):
    y = np.load(res_pairs_fnames[i])
    if not np.array_equal(x,y):
        print('Arrays not equal for i: %d' % i)

# first - distances that crossed the 4A threshold

# now we loop over every trajectory, and every distance
# we have 12720 distances total
traj0 = np.load(fnames[0])
dims = traj0.shape[1]
include_distance = [False] * dims

for fname in fnames:
    traj = np.load(fname)
    num_changed = (traj > 0.4).sum(0) * (traj < 0.4).sum(0)
    for i in range(dims):
        if num_changed[i] > 0:
            include_distance[i] = True

np.save('dist_cross/include_distance.npy', include_distance)

respairs_changed_indexes = []
for i in range(dims):
    if include_distance[i] == True:
        respairs_changed_indexes.append(i)

print(len(respairs_changed_indexes))
# 6567 (51.6%) distances classified

np.save('dist_cross/respairs_changed_indexes.npy', respairs_changed_indexes)

respairs = np.load(res_pairs_fnames[0])
respairs_changed = [respairs[i] for i in respairs_changed_indexes]

np.save('dist_cross/respairs_changed.npy', respairs_changed)

# now save trajectories with only these distances
x = glob.glob('data_cut_start_noH_stride10_featurized/dist/11707/*.npy')

for f in x:
    traj = np.load(f)[:,respairs_changed_indexes]
    np.save('data_cut_start_noH_stride10_featurized/dist_cross/11707/%s' % f[50:], traj)

x = glob.glob('data_cut_start_noH_stride10_featurized/dist/11709/*.npy')

for f in x:
    traj = np.load(f)[:,respairs_changed_indexes]
    np.save('data_cut_start_noH_stride10_featurized/dist_cross/11709/%s' % f[50:], traj)

# logistic (soft) distances
x = glob.glob('data_cut_start_noH_stride10_featurized/dist_soft/11707/*.npy')

for f in x:
    traj = np.load(f)[:,respairs_changed_indexes]
    np.save('data_cut_start_noH_stride10_featurized/dist_cross_soft/11707/%s' % f[55:], traj)
    
x = glob.glob('data_cut_start_noH_stride10_featurized/dist_soft/11709/*.npy')

for f in x:
    traj = np.load(f)[:,respairs_changed_indexes]
    np.save('data_cut_start_noH_stride10_featurized/dist_cross_soft/11709/%s' % f[55:], traj)    

# second - distances that crossed 4A and stayed for at least 10 ns
traj0 = np.load(fnames[0])
dims = traj0.shape[1]
include_distance_cross = [False] * dims
include_distance_time = [False] * dims

# cross
for fname in fnames:
    traj = np.load(fname)
    num_changed = (traj > 0.4).sum(0) * (traj < 0.4).sum(0)
    for i in range(dims):
        if num_changed[i] > 0:
            include_distance_cross[i] = True

# time
def find_distances_time(fname):
    include_distance_time_ = [False] * dims
    traj = np.load(fname)
    for i in range(dims):
        if include_distance_cross[i] == False:
            continue
        dim = traj[:,i]
        contact_formed = None
        for frame in dim:
            if frame < 0.4:
                if contact_formed is not None:
                    contact_formed += 1
                else:
                    contact_formed = 0
            else:
                if contact_formed is not None:
                    contact_formed = None
            if contact_formed == 1:
                include_distance_time_[i] = True
                break                     
    return include_distance_time_
    
pool = Pool(processes=18)
include_distance_time__ = pool.map(find_distances_time, fnames)
# save this intermediate list just in case
np.save('dist_time/include_distance_time__.npy', include_distance_time__)

# combine the info from all trajectories
for traj in include_distance_time__:
    for i in range(dims):
        if traj[i] == True:
            include_distance_time[i] = True

np.save('dist_time/include_distance_time.npy', include_distance_time)

respairs_changed_indexes_time = []
for i in range(dims):
    if include_distance_time[i] == True:
        respairs_changed_indexes_time.append(i)

print(len(respairs_changed_indexes_time))
# 5746 (45.1%) distances classified at 10 ns
# 5087 distances classified at 20 ns
# 4713 at 30 ns
# 4482 at 40 ns
# 4288 at 50 ns (33.7%)

np.save('dist_time/respairs_changed_indexes_time.npy', respairs_changed_indexes_time)

respairs = np.load(res_pairs_fnames[0])
respairs_changed_time = [respairs[i] for i in respairs_changed_indexes_time]

np.save('dist_time/respairs_changed_time.npy', respairs_changed_time)

# now save trajectories with only these distances
x = glob.glob('data_cut_start_noH_stride10_featurized/dist/11707/*.npy')

for f in x:
    traj = np.load(f)[:,respairs_changed_indexes_time]
    np.save('data_cut_start_noH_stride10_featurized/dist_time/11707/%s' % f[50:], traj)
    
x = glob.glob('data_cut_start_noH_stride10_featurized/dist/11709/*.npy')

for f in x:
    traj = np.load(f)[:,respairs_changed_indexes_time]
    np.save('data_cut_start_noH_stride10_featurized/dist_time/11709/%s' % f[50:], traj)

# logistic (soft) distances
x = glob.glob('data_cut_start_noH_stride10_featurized/dist_soft/11707/*.npy')

for f in x:
    traj = np.load(f)[:,respairs_changed_indexes_time]
    np.save('data_cut_start_noH_stride10_featurized/dist_time_soft/11707/%s' % f[55:], traj)
    
x = glob.glob('data_cut_start_noH_stride10_featurized/dist_soft/11709/*.npy')

for f in x:
    traj = np.load(f)[:,respairs_changed_indexes_time]
    np.save('data_cut_start_noH_stride10_featurized/dist_time_soft/11709/%s' % f[55:], traj)    
