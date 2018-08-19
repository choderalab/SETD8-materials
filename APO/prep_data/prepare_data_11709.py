import mdtraj as md
import glob

x = glob.glob('data/11709/*.h5')

i = 0
for f in x:
    print(i)
    i += 1
    traj = md.load(f)[50:]
    traj.save('data_cut_start/11709/%s' % f[11:])
    traj = traj.atom_slice(traj.top.select('not element H'))
    traj.save('data_cut_start_noH/11709/%s' % f[11:])
