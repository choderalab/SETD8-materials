import mdtraj as md
import os
import shutil

# put manual pdbs to be added in pdbs/
pdbs = ['4IJ8_inhibitor2.pdb', 'Apo_inhibitor2.pdb']

os.remove('templates/templates-resolved-seq.fa')
f = open('templates/templates-resolved-seq.fa', 'w')

for pdb in pdbs:
    shutil.copy(('pdbs/' + pdb), ('templates/structures-resolved/KMT5A_HUMAN_%s.pdb' % pdb.split('.')[0]))
    traj = md.load('pdbs/' + pdb)
    resolved_seq = traj.top.to_fasta()[0]
    f.write('\n>KMT5A_HUMAN_%s\n' % pdb.split('.')[0])
    f.write(resolved_seq)
    f.write('\n')
    
f.close()
    
    
