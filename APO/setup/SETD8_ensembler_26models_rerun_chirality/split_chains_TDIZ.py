import mdtraj as md

traj_a = md.load('files/TDIZ.pdb')
traj_b = md.load('files/TDIZ.pdb')

atoms_a = traj_a.top.select('chainid 0')
atoms_b = traj_a.top.select('chainid 1')

traj_a = traj_a.atom_slice(atoms_a)
traj_b = traj_b.atom_slice(atoms_b)

traj_a.save('manual_pdbs/TDIZ_A.pdb')
traj_b.save('manual_pdbs/TDIZ_B.pdb')
