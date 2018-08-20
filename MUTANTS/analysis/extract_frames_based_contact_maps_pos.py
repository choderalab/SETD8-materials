import numpy as np
import mdtraj as md
import multiprocessing
import glob

# list of contacts to check from extraction of top5
contacts_dict = {
 0: [(145, 83), (160, 108), (95, 56), (104, 80), (98, 37)],
 1: [(96, 79), (103, 24), (104, 54), (94, 71), (105, 80)],
 2: [(105, 80), (98, 37)],
 3: [(96, 79), (94, 55), (93, 55), (126, 55), (127, 55)],
 4: [(104, 96), (98, 37)],
 5: [(140, 35), (158, 35), (132, 47), (17, 14), (158, 107)],
 6: [(79, 75), (78, 74), (80, 76), (81, 78), (80, 77)],
 7: [(158, 68), (16, 13)],
 8: [(52, 25), (116, 54), (141, 51), (114, 51), (135, 51)],
 9: [(17, 14), (6, 3), (98, 37)],
 10: [(17, 14), (96, 79), (160, 77), (145, 116), (144, 116)],
 11: [(78, 71), (94, 71), (96, 76), (81, 78), (79, 76)],
 12: [],
 13: [(119, 8), (94, 71), (10, 7), (11, 7), (17, 14)],
 14: [(53, 23), (116, 54), (93, 55), (126, 55), (127, 55)],
 15: [(116, 54), (140, 35), (106, 102), (105, 100), (105, 101)],
 16: [(100, 29), (45, 23), (25, 22), (26, 23), (99, 55)],
 17: [(105, 80), (145, 109), (6, 3), (158, 107), (145, 106)],
 18: [(143, 79), (81, 78), (80, 77), (123, 120), (125, 120)],
 19: [(11, 8), (103, 24), (98, 37)],
 20: [(94, 71), (145, 80), (106, 37), (127, 84), (143, 80)],
 21: [(143, 54), (94, 55), (93, 55), (126, 55), (127, 55)],
 22: [(155, 83), (152, 83), (98, 37), (146, 116), (148, 145)],
 23: [(17, 14), (98, 37)],
 24: [(94, 71), (126, 55), (94, 81)],
 25: [(89, 60), (91, 60), (91, 59), (17, 14), (90, 60)],
 26: [(80, 76), (143, 79), (108, 79), (81, 78), (80, 77)]}

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
            if all([frame_contact_map[x[0],x[1]] for x in contacts_to_check]):
                frame_indexes.append((clone, frame_index))

    return frame_indexes

pool = multiprocessing.Pool(27)

frame_indexes = pool.map(get_run_contact_map, range(27))

np.save('extract_frames_indexes_all_top5_pos', frame_indexes)
