import numpy as np
import mdtraj as md
import multiprocessing
import glob

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

    return contact_maps

contact_maps = get_run_contact_map(21)

np.save('contact_maps_bool_run21', contact_maps)
