import pyemma
import mdtraj as md
import glob
import numpy as np

# data_cut_start_proteinonly_noH/11708/run0-clone211.h5 threw an error, so I'm excluding it
# also error on data_cut_start_proteinonly_noH/11710/run2-clone205.h5
# these files were very small - check later how many frames they have, if any?
fnames = [x for x in glob.glob('data_cut_start_proteinonly_noH/*/*.h5') if '11708/run0-clone211.h5' not in x and '11710/run2-clone205.h5' not in x]

traj = md.load(fnames[0])
top = traj.top

feat = pyemma.coordinates.featurizer(top)
feat.add_backbone_torsions(cossin = True)
feat.add_chi1_torsions(cossin = True)

source = pyemma.coordinates.source(fnames, features = feat)
X = source.get_output()

for i in range(len(X)):
    x = X[i]
    np.save('data_cut_start_proteinonly_noH_featurized/dih/%d.npy' % i, x)

np.save('data_cut_start_noH_stride10_featurized/dih_comb/X.npy', X)
