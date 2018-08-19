import mdtraj as md

# put manual pdbs to be added in manual_pdbs/
manual_pdbs = ['TDIX.pdb']

for pdb in manual_pdbs:
    traj = md.load('manual_pdbs/' + pdb)
    protein_atoms = traj.top.select('protein')
    traj = traj.atom_slice(protein_atoms)
    traj.save('templates/structures-resolved/SETD8_HUMAN_%s_A.pdb' % pdb.split('.')[0])
    resolved_seq = traj.top.to_fasta()[0]
    f = open('templates/templates-resolved-seq.fa', 'a')
    f.write('\n>SETD8_HUMAN_TDIX_A\n')
    f.write(resolved_seq)
    f.write('\n')
    f.close()
    
    
