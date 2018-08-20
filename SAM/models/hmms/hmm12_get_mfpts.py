import pyemma
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import glob

dtrajs = [np.load(x) for x in glob.glob('/cbio/jclab/home/rafal.wiewiora/repos/MSM_play/set8_ligands_11708_11710/dih/SAM/dtrajs/*.npy')]

hmm = pyemma.msm.estimate_hidden_markov_model(dtrajs, 12, 100, connectivity='largest', mincount_connectivity=25)

mfpts = np.zeros((hmm.nstates, hmm.nstates))

for i in range(hmm.nstates):
    for j in range(hmm.nstates):
        mfpts[i,j] = hmm.mfpt(i,j)

np.save('12/hmm12_mfpts', mfpts)        
