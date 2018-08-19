import glob
import numpy as np
import sklearn
import sklearn.model_selection

# we have kinetic / commute mapping - lag times: 1 and 10 - 4 sets of trajs

# let's first make splits by filename, then translate them to numbers for all 3: dihedrals, distances, log. distances

shuffle_split = sklearn.model_selection.ShuffleSplit(n_splits=10, test_size=0.5)
split_train = [[],[],[],[],[],[],[],[],[],[]]
split_test = [[],[],[],[],[],[],[],[],[],[]]

# 11707 has runs up to run168, 11709 up to run3
# 11707
for run in range(0,169):
    filenames_ = [x for x in glob.glob('/cbio/jclab/home/rafal.wiewiora/repos/MSM_play/set8_apo_11707_11709_FINAL/data_cut_start/11707/run%d-*.h5' % run) if 'run104-clone4' not in x]
    if filenames_:
        splits = np.array(list(shuffle_split.split(filenames_)))
        for index, split in enumerate(splits):
            for filename in split[0]:
                split_train[index].append(filenames_[filename].split('/')[-2] + '/' + filenames_[filename].split('/')[-1])
            for filename in split[1]:
                split_test[index].append(filenames_[filename].split('/')[-2] + '/' + filenames_[filename].split('/')[-1])

# 11709
for run in range(0,4):
    filenames_ = glob.glob('/cbio/jclab/home/rafal.wiewiora/repos/MSM_play/set8_apo_11707_11709_FINAL/data_cut_start/11709/run%d-*.h5' % run)
    if filenames_:
        splits = np.array(list(shuffle_split.split(filenames_)))
        for index, split in enumerate(splits):
            for filename in split[0]:
                split_train[index].append(filenames_[filename].split('/')[-2] + '/' + filenames_[filename].split('/')[-1])
            for filename in split[1]:
                split_test[index].append(filenames_[filename].split('/')[-2] + '/' + filenames_[filename].split('/')[-1])
        
np.save('splits_train_filenames', split_train)
np.save('splits_test_filenames', split_test)

# dihedrals

featurize_glob = glob.glob('/cbio/jclab/home/rafal.wiewiora/repos/MSM_play/set8_apo_11707_11709_FINAL/data_cut_start_noH/*/*.h5')
tica_glob = glob.glob('/cbio/jclab/home/rafal.wiewiora/repos/MSM_play/set8_apo_11707_11709_FINAL/data_cut_start_noH_stride10_featurized/dih/*.npy')

dihedrals_glob_dict = dict()
for i in range(len(tica_glob)):
    dihedrals_glob_dict[featurize_glob[int(tica_glob[i].split('/')[-1][:-4])].split('/')[-2] + '/' + featurize_glob[int(tica_glob[i].split('/')[-1][:-4])].split('/')[-1]] = i

split_train_dih = [[],[],[],[],[],[],[],[],[],[]]
split_test_dih = [[],[],[],[],[],[],[],[],[],[]]

for index, split in enumerate(split_train):
    for filename in split:
        split_train_dih[index].append(dihedrals_glob_dict[filename])

for index, split in enumerate(split_test):
    for filename in split:
        split_test_dih[index].append(dihedrals_glob_dict[filename])

np.save('splits_train_dih', split_train_dih)
np.save('splits_test_dih', split_test_dih)

# distances
tica_glob = glob.glob('/cbio/jclab/home/rafal.wiewiora/repos/MSM_play/set8_apo_11707_11709_FINAL/data_cut_start_noH_stride10_featurized/dist_cross/*/*.npy')

dist_glob_dict = dict()
for i in range(len(tica_glob)):
    dist_glob_dict[tica_glob[i].split('/')[-2] + '/' + tica_glob[i].split('/')[-1][:-3] + 'h5'] = i

split_train_dist = [[],[],[],[],[],[],[],[],[],[]]
split_test_dist = [[],[],[],[],[],[],[],[],[],[]]

for index, split in enumerate(split_train):
    for filename in split:
        split_train_dist[index].append(dist_glob_dict[filename])

for index, split in enumerate(split_test):
    for filename in split:
        split_test_dist[index].append(dist_glob_dict[filename])

np.save('splits_train_dist', split_train_dist)
np.save('splits_test_dist', split_test_dist)

# logistic distances
tica_glob = glob.glob('/cbio/jclab/home/rafal.wiewiora/repos/MSM_play/set8_apo_11707_11709_FINAL/data_cut_start_noH_stride10_featurized/dist_cross_soft/*/*.npy')

distlog_glob_dict = dict()
for i in range(len(tica_glob)):
    distlog_glob_dict[tica_glob[i].split('/')[-2] + '/' + tica_glob[i].split('/')[-1][:-3] + 'h5'] = i

split_train_distlog = [[],[],[],[],[],[],[],[],[],[]]
split_test_distlog = [[],[],[],[],[],[],[],[],[],[]]

for index, split in enumerate(split_train):
    for filename in split:
        split_train_distlog[index].append(distlog_glob_dict[filename])

for index, split in enumerate(split_test):
    for filename in split:
        split_test_distlog[index].append(distlog_glob_dict[filename])

np.save('splits_train_distlog', split_train_distlog)
np.save('splits_test_distlog', split_test_distlog)
