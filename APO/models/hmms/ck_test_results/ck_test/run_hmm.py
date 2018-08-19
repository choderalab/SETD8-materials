import pyemma
import numpy as np
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument('lag', type=int)
args = parser.parse_args()
lag = args.lag

dtrajs = list(np.load('/cbio/jclab/home/rafal.wiewiora/repos/MSM_play/set8_apo_11707_11709_FINAL/dih_full_dset/dtrajs_strided_fit/100.npy'))

hmm = pyemma.msm.estimate_hidden_markov_model(dtrajs, 25, lag, connectivity='largest')

np.save('%d/metastable_sets' % lag, hmm.metastable_sets)
np.save('%d/pi' % lag, hmm.pi)
np.save('%d/lifetimes' % lag, hmm.lifetimes)
np.save('%d/timescales' % lag, hmm.timescales())
np.save('%d/P' % lag, hmm.P)
np.save('%d/observation_probabilities' % lag, hmm.observation_probabilities)
np.save('%d/metastable_memberships' % lag, hmm.metastable_memberships)
