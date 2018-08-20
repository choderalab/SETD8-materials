import numpy as np
import mdtraj as md

indexes = np.load('mutant_12_corrected_new_strided_frame_indexes.npy')

traj = md.load('data/run12-clone%d.h5' % indexes[0][0])[indexes[0][1]]

for index in indexes[1:]:
    traj = md.join([traj, md.load('data/run12-clone%d.h5' % index[0])[index[1]]])
    
traj.save('mutant_12_corrected_new_frames_traj.h5')    
