import pyemma
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

dtrajs = list(np.load('/cbio/jclab/home/rafal.wiewiora/repos/MSM_play/set8_apo_11707_11709_FINAL/dih_full_dset/dtrajs_strided_fit/100.npy'))

hmm = pyemma.msm.estimate_hidden_markov_model(dtrajs, 25, 100, connectivity='largest')

mfpts = np.zeros((hmm.nstates, hmm.nstates))

for i in range(hmm.nstates):
    for j in range(hmm.nstates):
        mfpts[i,j] = hmm.mfpt(i,j)
        
np.save('hmms/hmm25_mfpts', mfpts)        
    
