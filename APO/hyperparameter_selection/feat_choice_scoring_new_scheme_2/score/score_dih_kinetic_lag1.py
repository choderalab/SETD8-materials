import numpy as np
import pyemma

lag = 1
mapping = 'kinetic'
feat = 'dih'

cluster_nos = [50, 100, 500, 1000]
replicates = [0,1,2,3,4]

pyemma_scores_1 = []
pyemma_scores_1_train = []

pyemma_scores_2 = []
pyemma_scores_2_train = []

for cluster_no in cluster_nos:
    pyemma_scores_1.append([])
    pyemma_scores_1_train.append([])
    
    pyemma_scores_2.append([])
    pyemma_scores_2_train.append([])
    
    for replicate in replicates:
        dtrajs_train = list(np.load('/data/chodera/rafal.wiewiora/set8_apo_11707_11709/%s/tica_trajs_%s_lag%d/dtrajs/%d_%d_train.npy' % (feat, mapping, lag, cluster_no, replicate)))
        dtrajs_test = list(np.load('/data/chodera/rafal.wiewiora/set8_apo_11707_11709/%s/tica_trajs_%s_lag%d/dtrajs/%d_%d_test.npy' % (feat, mapping, lag, cluster_no, replicate)))
        dtrajs_test = [np.concatenate(x) for x in dtrajs_test]

        try:
            # VAMP-1
            pyemma_msm = pyemma.msm.estimate_markov_model(dtrajs_train, lag=10, score_method='VAMP1', score_k=10)
            pyemma_score_1 = pyemma_msm.score(dtrajs_test, score_method='VAMP1', score_k=10)
            pyemma_score_1_train = pyemma_msm.score(dtrajs_train, score_method='VAMP1', score_k=10)
        except:
            pyemma_score_1 = 'error'
            pyemma_score_1_train = 'error'
            
        try:    
            # VAMP-2
            pyemma_msm = pyemma.msm.estimate_markov_model(dtrajs_train, lag=10, score_method='VAMP2', score_k=10)
            pyemma_score_2 = pyemma_msm.score(dtrajs_test, score_method='VAMP2', score_k=10)
            pyemma_score_2_train = pyemma_msm.score(dtrajs_train, score_method='VAMP2', score_k=10)
        except:
            pyemma_score_2 = 'error'
            pyemma_score_2_train = 'error'
        
        pyemma_scores_1[-1].append(pyemma_score_1)
        pyemma_scores_1_train[-1].append(pyemma_score_1_train)
        
        pyemma_scores_2[-1].append(pyemma_score_2)
        pyemma_scores_2_train[-1].append(pyemma_score_2_train)
        
        print((pyemma_score_1, pyemma_score_1_train, pyemma_score_2, pyemma_score_2_train))

pyemma_scores = [pyemma_scores_1, pyemma_scores_1_train, pyemma_scores_2, pyemma_scores_2_train]

np.save('pyemma_scores_%s_%s_lag%d' % (feat, mapping, lag), pyemma_scores) 
