import numpy as np
import glob
import mdtraj as md

wt_combined_traj = md.load('/cbio/jclab/home/rafal.wiewiora/repos/MSM_play/set8_apo_cancer_mutants_11713/all_wt_CA_only_traj.h5')  

for run in glob.glob('11713_extracted_frames_pos_newconf/*.h5'):
    traj = md.load(run)
    rmsds = []
    for frame in range(len(traj)):
        rmsd = md.rmsd(wt_combined_traj, traj, frame=frame)
        rmsds.append(rmsd)
    np.save('11713_extracted_frames_pos_newconf_rmsds/' + run.split('/')[1][:-3], rmsds)        

for run in glob.glob('11713_extracted_frames_neg_newconf/*.h5'):
    # run12 separately because have to remove the missing CA from wt traj
    if 'run12' in run:
        continue
    traj = md.load(run)
    rmsds = []
    for frame in range(len(traj)):
        rmsd = md.rmsd(wt_combined_traj, traj, frame=frame)
        rmsds.append(rmsd)
    np.save('11713_extracted_frames_neg_newconf_rmsds/' + run.split('/')[1][:-3], rmsds)
    
# run 12 only - modify the WT - remove residue 92
wt_combined_traj = wt_combined_traj.atom_slice(wt_combined_traj.top.select('not resid 92'))

for run in glob.glob('11713_extracted_frames_neg_newconf/*.h5'):
    # run12 separately because have to remove the missing CA from wt traj
    if not ('run12' in run):
        continue
    traj = md.load(run)
    rmsds = []
    for frame in range(len(traj)):
        rmsd = md.rmsd(wt_combined_traj, traj, frame=frame)
        rmsds.append(rmsd)
    np.save('11713_extracted_frames_neg_newconf_rmsds/' + run.split('/')[1][:-3], rmsds)

# extract those frames from the dataset, save as h5, dcd and pdb
# create a dictionary of starting indexes in the wt_traj vs trajectory path
wt_indexes_dict = dict()
wt_trajs = glob.glob('/cbio/jclab/home/rafal.wiewiora/repos/MSM_play/set8_apo_11707_11709_FINAL/data_cut_start_noH_stride10/*/*.h5')

# the wt_combined_traj has the first traj twice by mistake - repeat that here
# lengths are one index bigger than the last index of the given trajectory - so we look for the traj that first
# has the index higher than the frame we're looking for and that's the one
current_length = len(md.open(wt_trajs[0]))
wt_indexes_dict = dict()
lengths = []
lengths.append(current_length)
for traj in wt_trajs:
    current_length += len(md.open(traj))
    wt_indexes_dict[current_length] = traj
    lengths.append(current_length)
    print(current_length)
    
lengths = np.sort(lengths) 

# extract frames
# pos
for run in glob.glob('11713_extracted_frames_pos_newconf_rmsds/*.npy'):
    run_ = np.load(run)
    for i, frame in enumerate(run_):
        # find the trajectory
        for j, length in enumerate(lengths):
            if length > np.argmin(frame):
                traj = wt_indexes_dict[length]
                frame_index = np.argmin(frame) - lengths[j-1]
                break
        traj = md.load(traj)[frame_index]
        if i == 0:
            traj_combined = traj
        else:
            traj_combined = md.join([traj_combined, traj])
    traj_combined.save('11713_extracted_frames_pos_newconf_rmsds/%s.h5' % run.split('/')[-1][:-4])        
    traj_combined.save('11713_extracted_frames_pos_newconf_rmsds/%s.dcd' % run.split('/')[-1][:-4])
    
# neg
for run in glob.glob('11713_extracted_frames_neg_newconf_rmsds/*.npy'):
    run_ = np.load(run)
    for i, frame in enumerate(run_):
        # find the trajectory
        for j, length in enumerate(lengths):
            if length > np.argmin(frame):
                traj = wt_indexes_dict[length]
                frame_index = np.argmin(frame) - lengths[j-1]
                break
        traj = md.load(traj)[frame_index]
        if i == 0:
            traj_combined = traj
        else:
            traj_combined = md.join([traj_combined, traj])
    traj_combined.save('11713_extracted_frames_neg_newconf_rmsds/%s.h5' % run.split('/')[-1][:-4])
    traj_combined.save('11713_extracted_frames_neg_newconf_rmsds/%s.dcd' % run.split('/')[-1][:-4])    
