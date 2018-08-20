import mdtraj as md
import glob
from multiprocessing import Pool

x = glob.glob('data/11708/*.h5')
n_threads = 10

def prepare_data(f):
    traj = md.load(f)[50:]
    traj.save('data_cut_start/11708/%s' % f.split('/')[-1])
    traj = traj.atom_slice(traj.top.select('protein and not element H'))
    if len(list(traj.top.chains)) > 1 and len(list(traj.top.chains)) == 2:
        traj = traj.atom_slice(traj.top.select('chainid 0'))
    elif len(list(traj.top.chains)) > 1 and len(list(traj.top.chains)) != 2:
        raise Exception('Too many chains in %s' % f)
    traj.save('data_cut_start_proteinonly_noH/11708/%s' % f.split('/')[-1])

pool = Pool(n_threads)
pool.map(prepare_data, x)
