import mdtraj as md

manual_pdbs = ['TDIZ_A.pdb', 'TDIZ_B.pdb']

f = open('templates/templates-resolved-seq.fa', 'w')

for pdb in manual_pdbs:
    traj = md.load('manual_pdbs/' + pdb)
    protein_atoms = traj.top.select('protein')
    traj = traj.atom_slice(protein_atoms)
    traj.save('templates/structures-resolved/KMT5A_HUMAN_%s.pdb' % pdb.split('.')[0])
    resolved_seq = traj.top.to_fasta()[0]
    f.write('\n>KMT5A_HUMAN_%s\n' % pdb.split('.')[0])
    f.write(resolved_seq)
    f.write('\n')
    
f.close()
