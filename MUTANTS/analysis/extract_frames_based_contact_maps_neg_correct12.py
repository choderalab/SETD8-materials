import numpy as np
import mdtraj as md
import multiprocessing
import glob

# list of contacts to check from extraction of top5
# subtracting one from any index over 92
contacts_dict = {
 12: [(92, 81), (81, 92), (93, 80), (80, 93), (80, 94)],
}

# now with Booleyan contacts
def get_run_contact_map(run):

    frame_indexes = []

    for clone in range(40):
        traj = md.load('data/run%d-clone%d.h5' % (run,clone))
        # trajs are 0.5 ns / frame, we want only the part after 750 ns == 1500 frames
        if not (len(traj) > 1500):
            continue
        else:
            traj = traj[1500:]
        distances, residue_pairs = md.compute_contacts(traj)
        contact_map = md.geometry.squareform(distances, residue_pairs)
        contact_map_bool = contact_map < 0.4

        contacts_to_check = contacts_dict[run]
        for frame_index, frame_contact_map in enumerate(contact_map_bool):
            if not any([frame_contact_map[x[0],x[1]] for x in contacts_to_check]):
                frame_indexes.append((clone, frame_index))

    return frame_indexes

#pool = multiprocessing.Pool(27)

frame_indexes = get_run_contact_map(12)

np.save('extract_frames_indexes_all_top5_neg_correct12', frame_indexes)
