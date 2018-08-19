import pyemma
import numpy as np
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument('no_macrostates', type=int)
args = parser.parse_args()
no_macrostates = args.no_macrostates

dtrajs = list(np.load('100.npy'))

hmm = pyemma.msm.estimate_hidden_markov_model(dtrajs, no_macrostates, 100, connectivity='largest')

np.save('%d/metastable_sets' % no_macrostates, hmm.metastable_sets)
np.save('%d/pi' % no_macrostates, hmm.pi)
np.save('%d/lifetimes' % no_macrostates, hmm.lifetimes)
np.save('%d/timescales' % no_macrostates, hmm.timescales())
np.save('%d/P' % no_macrostates, hmm.P)
np.save('%d/observation_probabilities' % no_macrostates, hmm.observation_probabilities)
np.save('%d/metastable_memberships' % no_macrostates, hmm.metastable_memberships)
