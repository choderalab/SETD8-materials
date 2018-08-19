import pyemma
import mdtraj as md
import glob
import numpy as np

fnames = glob.glob('data_cut_start_noH_stride10/*/*.h5')

traj = md.load(fnames[0])
top = traj.top

feat = pyemma.coordinates.featurizer(top)
feat.add_backbone_torsions(cossin = True)
feat.add_chi1_torsions(cossin = True)

source = pyemma.coordinates.source(fnames, features = feat)
X = source.get_output()

for i in range(len(X)):
    x = X[i]
    np.save('data_cut_start_noH_stride10_featurized/dih/%d.npy' % x, x)

np.save('data_cut_start_noH_stride10_featurized/dih_comb/X.npy', X)
