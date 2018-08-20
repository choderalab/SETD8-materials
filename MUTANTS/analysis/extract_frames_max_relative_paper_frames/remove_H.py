import mdtraj as md

for run in [3,5,6,10,11,12,13,14,15,17,18,21,22]:
    traj = md.load('cluster_traj_%d.pdb' % run)
    traj = traj.atom_slice(traj.top.select('not element H'))
    traj.save('cluster_traj_%d_noH.pdb' % run)
