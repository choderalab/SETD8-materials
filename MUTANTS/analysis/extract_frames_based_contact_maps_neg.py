import numpy as np
import mdtraj as md
import multiprocessing
import glob

# list of contacts to check from extraction of top5
contacts_dict = {
 0: [(7, 4), (160, 153), (105, 101), (160, 157), (156, 152)],
 1: [(16, 12), (120, 12), (15, 12), (119, 12), (159, 145)],
 2: [(158, 80), (159, 145), (156, 83), (158, 79), (159, 153)],
 3: [(93, 56), (101, 55), (93, 54), (95, 55), (104, 55)],
 4: [(159, 145), (159, 108), (158, 143), (42, 27), (159, 153)],
 5: [(158, 79), (159, 145), (156, 81), (158, 81), (156, 152)],
 6: [(79, 71), (158, 80), (96, 79), (20, 1), (16, 1)],
 7: [(49, 45), (158, 80), (50, 19), (17, 12), (156, 81)],
 8: [(106, 52), (50, 45), (159, 108), (159, 145), (160, 153)],
 9: [(65, 61), (83, 65), (66, 62), (83, 68), (83, 64)],
 10: [(159, 145), (159, 108), (80, 54), (158, 79), (81, 54)],
 11: [(78, 74), (158, 79), (79, 71), (158, 81), (156, 81)],
 12: [(93, 81), (81, 93), (94, 80), (80, 94), (80, 95)],
 13: [(148, 144), (152, 147), (98, 37), (104, 98), (159, 145)],
 14: [(101, 24), (101, 53), (102, 24), (101, 22), (93, 56)],
 15: [(101, 24), (101, 53), (101, 22), (102, 24), (52, 24)],
 16: [(101, 22), (119, 16), (53, 19), (103, 24), (126, 19)],
 17: [(159, 145), (144, 114), (144, 116), (144, 115), (156, 83)],
 18: [(123, 86), (79, 71), (120, 12), (118, 85), (158, 80)],
 19: [(160, 157), (15, 12), (80, 37), (104, 96), (16, 12)],
 20: [(158, 81), (158, 79), (159, 145), (158, 143), (159, 108)],
 21: [(101, 55), (93, 56), (93, 54), (158, 79), (116, 82)],
 22: [(159, 108), (158, 108), (119, 16), (158, 80), (158, 79)],
 23: [(158, 80), (119, 16), (156, 81), (160, 157), (156, 152)],
 24: [(89, 60), (158, 79), (159, 145), (148, 110), (80, 37)],
 25: [(16, 12), (15, 12), (49, 14), (17, 13), (120, 12)],
 26: [(159, 108), (158, 80), (79, 71), (159, 145), (96, 79)]}

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

pool = multiprocessing.Pool(27)

frame_indexes = pool.map(get_run_contact_map, range(27))

np.save('extract_frames_indexes_all_top5_neg', frame_indexes)
