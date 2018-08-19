import mdtraj as md
import glob

x = glob.glob('data_cut_start_noH/11709/*.h5')

i = 0
for f in x:
    print(i)
    i += 1
    traj = md.load(f)[::10]
    traj.save('data_cut_start_noH_stride10/11709/%s' % f[25:])
