import mdtraj as md
import os

traj_a = md.load('files/TDIY.pdb')
traj_b = md.load('files/TDIY.pdb')

atoms_a = traj_a.top.select('chainid 0')
atoms_b = traj_a.top.select('chainid 1')

traj_a.restrict_atoms(atoms_a)
traj_b.restrict_atoms(atoms_b)

os.mkdir('manual_pdbs')
traj_a.save('manual_pdbs/TDIY_A.pdb')
traj_b.save('manual_pdbs/TDIY_B.pdb')
