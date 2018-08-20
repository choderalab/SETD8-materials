import numpy as np
import mdtraj as md
import multiprocessing
import glob
import pyemma

# list of contacts to check from extraction of top5
contacts_dict = {11: [(76, 71),
  (76, 79),
  (76, 96),
  (158, 35),
  (158, 107),
  (75, 156),
  (78, 81),
  (106, 145)]}
 
#frame_indexes = np.load('extract_frames_max_relative_paper_run11.npy')

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
    #frames = frame_indexes
    # now cluster
    all_frames_traj = md.load('all_frames_traj_%d.h5' % run)
    top = all_frames_traj.top
    #del all_frames_traj
    data = pyemma.coordinates.load('all_frames_traj_%d.h5' % run, top=top)
    clustering = pyemma.coordinates.cluster_regspace(data, dmin=0.3, metric='minRMSD')
    cluster_indexes = clustering.index_clusters

    cluster_indexes_ = []
    for cluster in cluster_indexes:
        center = cluster[0][1]
        cluster_indexes_.append(center)
        
    frame = cluster_indexes_[0]
    cluster_traj = all_frames_traj[frame]
    for frame in cluster_indexes_[1:]:
        cluster_traj = md.join([cluster_traj,all_frames_traj[frame]])
    
    return cluster_traj

# now run stuff
#all_frames_traj = get_run_contact_map(11)

#all_frames_traj.save('all_frames_traj_11.h5')

cluster_traj = cluster_trajs(11)

cluster_traj.save('cluster_traj_11_approach2.h5')
cluster_traj.save('cluster_traj_11_approach2.pdb')
