import mdtraj as md
import glob

x = glob.glob('data/11707/*.h5')

broken_runs = [1,7,8,10,16,22,24,25,26,27,28,30,31,32,33,34,35,39,40,41,60,61,62,63,64,65,66,67,68,69,75,76,77,78,79]

i = 0
for f in x:
    print(i)
    i += 1
    if not (int(f.split('-')[0][14:])) in broken_runs:
        traj = md.load(f)[50:]
        traj.save('data_cut_start/11707/%s' % f[11:])
        traj = traj.atom_slice(traj.top.select('not element H'))
        traj.save('data_cut_start_noH/11707/%s' % f[11:])
