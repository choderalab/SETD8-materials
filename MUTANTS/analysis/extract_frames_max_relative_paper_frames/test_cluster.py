import numpy as np
import mdtraj as md
import multiprocessing
import glob
import pyemma

# list of contacts to check from extraction of top5
contacts_dict = {3: [(55, 93), (55, 94), (55, 126), (55, 127)],
 5: [(134, 1), (158, 35), (158, 107)],
 6: [(134, 1)],
 10: [(152, 82), (157, 80), (160, 74), (160, 77)],
 11: [(76, 71),
   (76, 79),
   (76, 96),
   (158, 35),
   (158, 107),
   (75, 156),
   (78, 81),
   (106, 145)],
 12: [(82, 93),
  (89, 58),
  (89, 59),
  (90, 58),
  (90, 92),
  (91, 57),
  (91, 93),
  (83, 92),
  (61, 85),
  (61, 88)],
 13: [(8, 119), (15, 4), (48, 3), (48, 5), (68, 149), (72, 151), (79, 159)],
 14: [(55, 58), (55, 100), (100, 97), (101, 54), (102, 54)],
 15: [(100, 37)],
 17: [(106, 145), (158, 107)],
 18: [(124, 119)],
 21: [(55, 58), (55, 93), (55, 94), (55, 126), (55, 127), (158, 107)],
 22: [(146, 116), (146, 82)]}
 
# test that all frames for one of the runs have the required contacts established

def test_run(run):
    traj = md.load('cluster_traj_%d.pdb' % run)
    distances, residue_pairs = md.compute_contacts(traj)
    contact_map = md.geometry.squareform(distances, residue_pairs)
    contact_map_bool = contact_map < 0.4

    contacts_to_check = contacts_dict[run]
    for frame_index, frame_contact_map in enumerate(contact_map_bool):
        if not all([frame_contact_map[x[0],x[1]] for x in contacts_to_check]):
            raise Exception('Test fail: run %d, frame %d' % (run, frame_index))
            
    print('Run %d pass!' % run)            
            
for run in contacts_dict:
    test_run(run)            
