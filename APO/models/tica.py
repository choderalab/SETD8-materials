import pyemma
import glob
import numpy as np
import mdtraj as md

lag_time = 1

# remove the unfolding trajectory - number 3959
source = pyemma.coordinates.source([x for x in glob.glob('../data_cut_start_noH_stride10_featurized/dih/*.npy') if '3959' not in x])

tica_commute = pyemma.coordinates.tica(lag=lag_time, kinetic_map=False, commute_map=True, var_cutoff=1)
stages = [source, tica_commute]
pipeline = pyemma.coordinates.pipeline(stages, chunksize = 1000)

eigenvalues = tica_commute.eigenvalues
eigenvectors = tica_commute.eigenvectors
timescales = tica_commute.timescales
correlation = tica_commute.feature_TIC_correlation

np.save('eigenvalues', eigenvalues)
np.save('eigenvectors', eigenvectors)
np.save('timescales', timescales)
np.save('correlation', correlation)

trajs = [x for x in glob.glob('../data_cut_start_noH_featurized/dih/*.npy') if '3959' not in x]

for traj in trajs:
    name = traj.split('/')[3]
    traj_ = tica_commute.transform(np.load(traj))
    np.save('trajs/%s' % name, traj_)

# we saved all dimensions, let's cut trajs down to 95% kinetic content - up to including 466
tica_trajs = glob.glob('trajs/*.npy')

for tica_traj in tica_trajs:
    name = tica_traj.split('/')[1]
    traj_ = np.load(tica_traj)
    traj_ = traj_[:, :467]
    np.save('trajs_95/%s' % name, traj_)
    
# also save just top 10 TICS for analysis
for tica_traj in tica_trajs:
    name = tica_traj.split('/')[1]
    traj_ = np.load(tica_traj)
    traj_ = traj_[:, :10]
    np.save('trajs_10/%s' % name, traj_)

# transform starting structures
all_traj = md.load('starting_structures/all_traj.h5')
top = all_traj.top

feat = pyemma.coordinates.featurizer(top)
feat.add_backbone_torsions(cossin = True)
feat.add_chi1_torsions(cossin = True)

source = pyemma.coordinates.source('starting_structures/all_traj.h5', features = feat)
X = source.get_output()
Y = tica_commute.transform(X)
np.save('starting_structures/starting_structures.npy', Y)
