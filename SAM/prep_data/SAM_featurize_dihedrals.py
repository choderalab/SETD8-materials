import pyemma
import mdtraj as md
import glob
import numpy as np

fnames = glob.glob('data_cut_start_proteinonly_noH/11710/run0-*.h5')  + [x for x in glob.glob('data_cut_start_proteinonly_noH/11710/run2-*.h5') if '11710/run2-clone205.h5' not in x]

traj = md.load(fnames[0])
top = traj.top

# check that all 999 trajectories have the same topology
for fname in fnames:
    traj_ = md.load(fname)
    top_ = traj_.top
    if top_ != top:
        print(fname)

feat = pyemma.coordinates.featurizer(top)
feat.add_backbone_torsions(cossin = True)
feat.add_chi1_torsions(cossin = True)

source = pyemma.coordinates.source(fnames, features = feat)
X = source.get_output()

for i in range(len(X)):
    x = X[i]
    np.save('data_cut_start_proteinonly_noH_featurized/dih/SAM/%d.npy' % i, x)
