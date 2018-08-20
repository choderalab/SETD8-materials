import numpy as np
import mdtraj as md
import multiprocessing
import glob

# list of contacts to check from extraction of top5
contacts_dict = {3: [(55, 93), (55, 94), (55, 126), (55, 127)],
 5: [(134, 1), (158, 35), (158, 107)],
 6: [(134, 1)],
 8: [(76, 71),
  (76, 79),
  (76, 96),
  (158, 35),
  (158, 107),
  (75, 156),
  (78, 81),
  (106, 145)],
 10: [(152, 82), (157, 80), (160, 74), (160, 77)],
 12: [(80, 71),
  (82, 93),
  (82, 94),
  (82, 157),
  (89, 58),
  (89, 59),
  (92, 55),
  (90, 58),
  (90, 59),
  (90, 64),
  (90, 83),
  (90, 92),
  (91, 57),
  (91, 93),
  (91, 124),
  (91, 123),
  (83, 92),
  (83, 67),
  (95, 67),
  (61, 85),
  (61, 88)],
 13: [(8, 119), (15, 4), (48, 3), (48, 5), (68, 149), (72, 151), (79, 159)],
 14: [(55, 58),
  (55, 93),
  (55, 94),
  (55, 100),
  (55, 126),
  (55, 127),
  (23, 19),
  (23, 52),
  (23, 53),
  (100, 97),
  (101, 54),
  (101, 95),
  (101, 97),
  (102, 54)],
 15: [(100, 37)],
 17: [(106, 145), (158, 107)],
 18: [(124, 119)],
 21: [(55, 58), (55, 93), (55, 94), (55, 126), (55, 127), (158, 107)],
 22: [(146, 116), (146, 82)]}

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

pool = multiprocessing.Pool(13)

frame_indexes = pool.map(get_run_contact_map, [3,5,6,8,10,12,13,14,15,17,18,21,22])

np.save('extract_frames_max_relative_paper', frame_indexes)
