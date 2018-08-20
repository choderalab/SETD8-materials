import numpy as np
import mdtraj as md
import multiprocessing
import glob

def get_run_contact_map(run):

    contact_maps = []

    for clone in range(40):
        traj = md.load('data/run%d-clone%d.h5' % (run, clone))
        # trajs are 0.5 ns / frame, we want only the part after 750 ns == 1500 frames
        if not (len(traj) > 1500):
            continue
        else:
            traj = traj[1500:]
        distances, residue_pairs = md.compute_contacts(traj)
        contact_map = md.geometry.squareform(distances, residue_pairs)
        contact_map = np.mean(contact_map, axis=0)
        contact_maps.append(contact_map)

    contact_map = np.mean(np.array(contact_maps), axis=0)

    return contact_map

pool = multiprocessing.Pool(32)

contact_maps = pool.map(get_run_contact_map, range(27))

np.save('contact_maps', contact_maps)

# now with Booleyan contacts
def get_run_contact_map(run):

    contact_maps = []

    for clone in range(40):
        traj = md.load('data/run%d-clone%d.h5' % (run, clone))
        # trajs are 0.5 ns / frame, we want only the part after 750 ns == 1500 frames
        if not (len(traj) > 1500):
            continue
        else:
            traj = traj[1500:]
        distances, residue_pairs = md.compute_contacts(traj)
        contact_map = md.geometry.squareform(distances, residue_pairs)
        contact_map_bool = contact_map < 0.4
        contact_map_bool_float = contact_map_bool.astype('float')
        contact_map = np.mean(contact_map_bool_float, axis=0)
        contact_maps.append(contact_map)

    contact_map = np.mean(np.array(contact_maps), axis=0)

    return contact_map

pool = multiprocessing.Pool(32)

contact_maps = pool.map(get_run_contact_map, range(27))

np.save('contact_maps_bool', contact_maps)

# now get the same for the WT data started from 1ZKK_A (p11707, run0 and run133) as control
trajs = glob.glob('/cbio/jclab/home/rafal.wiewiora/repos/MSM_play/set8_apo_11707_11709_FINAL/data/11707/run0-*') + glob.glob('/cbio/jclab/home/rafal.wiewiora/repos/MSM_play/set8_apo_11707_11709_FINAL/data/11707/run133-*')

contact_maps = []

for traj_ in trajs:                                           
    traj = md.load(traj_)
    # trajs are 0.5 ns / frame, we want only the part after 750 ns == 1500 frames
    if not (len(traj) > 1500):
        continue
    else:
        traj = traj[1500:]
    distances, residue_pairs = md.compute_contacts(traj)
    contact_map = md.geometry.squareform(distances, residue_pairs)
    contact_map = np.mean(contact_map, axis=0)
    contact_maps.append(contact_map)

contact_map = np.mean(np.array(contact_maps), axis=0)

np.save('wt_contact_map', contact_map)

# now with Boolean contacts
contact_maps = []

for traj_ in trajs:
    traj = md.load(traj_)
    # trajs are 0.5 ns / frame, we want only the part after 750 ns == 1500 frames
    if not (len(traj) > 1500):
        continue
    else:
        traj = traj[1500:]
    distances, residue_pairs = md.compute_contacts(traj)
    contact_map = md.geometry.squareform(distances, residue_pairs)
    contact_map_bool = contact_map < 0.4
    contact_map_bool_float = contact_map_bool.astype('float')
    contact_map = np.mean(contact_map_bool_float, axis=0)
    contact_maps.append(contact_map)

contact_map = np.mean(np.array(contact_maps), axis=0)

np.save('wt_contact_map_bool', contact_map)
