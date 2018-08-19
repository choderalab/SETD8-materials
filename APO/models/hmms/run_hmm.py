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

dtrajs = list(np.load('/cbio/jclab/home/rafal.wiewiora/repos/MSM_play/set8_apo_11707_11709_FINAL/dih_full_dset/dtrajs_strided_fit/100.npy'))

hmm = pyemma.msm.estimate_hidden_markov_model(dtrajs, no_macrostates, 100, connectivity='largest')

np.save('%d/metastable_sets' % no_macrostates, hmm.metastable_sets)
np.save('%d/pi' % no_macrostates, hmm.pi)
np.save('%d/lifetimes' % no_macrostates, hmm.lifetimes)
np.save('%d/timescales' % no_macrostates, hmm.timescales())
np.save('%d/P' % no_macrostates, hmm.P)
np.save('%d/observation_probabilities' % no_macrostates, hmm.observation_probabilities)
np.save('%d/metastable_memberships' % no_macrostates, hmm.metastable_memberships)

for i in range(len(hmm.observation_probabilities)):
    plt.figure(dpi=300)
    plt.scatter(range(len(hmm.observation_probabilities[i])), hmm.observation_probabilities[i])
    plt.plot(range(len(hmm.observation_probabilities[i])), hmm.observation_probabilities[i])
    plt.savefig('%d/observation_probabilities/%s.png' % (no_macrostates, i))
    plt.close()
    
for i in range(len(hmm.metastable_memberships)):
    plt.figure(dpi=300)
    plt.plot(range(len(hmm.metastable_memberships[i])), hmm.metastable_memberships[i])
    plt.scatter(range(len(hmm.metastable_memberships[i])), hmm.metastable_memberships[i])
    plt.savefig('%d/metastable_memberships/%s.png' % (no_macrostates,i))
    plt.close()
