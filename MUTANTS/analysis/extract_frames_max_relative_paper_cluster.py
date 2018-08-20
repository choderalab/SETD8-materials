import numpy as np
import mdtraj as md
import multiprocessing
import glob
import pyemma

# list of contacts to check from extraction of top5
contacts_dict = {3: [(55, 93), (55, 94), (55, 126), (55, 127)],
 5: [(134, 1), (158, 35), (158, 107)],
 6: [(134, 1)],
 8: [(76, 96)],
 10: [(152, 82), (157, 80), (160, 74), (160, 77)],
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
 
frame_indexes = np.load('extract_frames_max_relative_paper_fixnoframes.npy')

# now with Booleyan contacts
def get_run_contact_map(run):
    
    traj_started = False

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
                if traj_started:
                    all_frames_traj = md.join([all_frames_traj, traj[frame_index].atom_slice(traj.top.select('name CA'))])
                else:
                    all_frames_traj = traj[frame_index].atom_slice(traj.top.select('name CA'))
                    traj_started = True
                    
    return all_frames_traj  
    
def cluster_trajs(run):
    no = {3:0,5:1,6:2,8:3,10:4,12:5,13:6,14:7,15:8,17:9,18:10,21:11,22:12}[run]
    frames = frame_indexes[no]
    # now cluster
    all_frames_traj = md.load('all_frames_traj_%d.h5' % run)       
    top = all_frames_traj.top
    del all_frames_traj
    data = pyemma.coordinates.load('all_frames_traj_%d.h5' % run, top=top)
    clustering = pyemma.coordinates.cluster_regspace(data, dmin=0.3, metric='minRMSD')
    cluster_indexes = clustering.index_clusters
    
    cluster_indexes_ = []
    for cluster in cluster_indexes:
        center = cluster[0][1]
        cluster_indexes_.append(frames[center])
        
    frame = cluster_indexes_[0]
    cluster_traj = md.load('data/run%d-clone%d.h5' % (run,frame[0]))[frame[1]+1500]
    for frame in cluster_indexes_[1:]:
        cluster_traj = md.join([cluster_traj,md.load('data/run%d-clone%d.h5' % (run,frame[0]))[frame[1]+1500]])

    return cluster_traj

# now run stuff
pool = multiprocessing.Pool(13)

all_frames_trajs = pool.map(get_run_contact_map, [3,5,6,8,10,12,13,14,15,17,18,21,22])

for i,traj in enumerate(all_frames_trajs):
    no = [3,5,6,8,10,12,13,14,15,17,18,21,22][i]
    traj.save('all_frames_traj_%d.h5' % no)
    
pool = multiprocessing.Pool(13) 
    
cluster_trajs = pool.map(cluster_trajs, [3,5,6,8,10,12,13,14,15,17,18,21,22])

for i,traj in enumerate(cluster_trajs):
    no = [3,5,6,8,10,12,13,14,15,17,18,21,22][i]
    traj.save('cluster_traj_%d.h5' % no)
    traj.save('cluster_traj_%d.pdb' % no)
