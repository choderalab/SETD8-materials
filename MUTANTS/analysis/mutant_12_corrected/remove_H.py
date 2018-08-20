import mdtraj as md

traj = md.load('mutant_12_corrected_new_frames_traj.pdb')
traj = traj.atom_slice(traj.top.select('not element H'))
traj.save('mutant_12_corrected_new_frames_traj_noH.pdb')
