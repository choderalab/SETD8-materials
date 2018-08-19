import mdtraj as md
import os
# put manual pdbs to be added in manual_pdbs/
manual_pdbs = ['TDIY_A.pdb', 'TDIY_B.pdb']

# special to this version - we couldn't add the 24th and 25th models to the ready 23 set (clustering doesn't recognize new addit# ions) have to delete all template sequences and only paste in the extra one we want to do here
os.remove('templates/templates-resolved-seq.fa')
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
    
    
