import mdtraj as md
import os
import shutil
# put manual pdbs to be added in pdbs/
pdbs = ['1ZKK_apo.pdb', '4IJ8_apo.pdb', 'Inhibitor_apo.pdb', 'Inhibitor_4IJ8.pdb']

# special to this version - we couldn't add the 24th and 25th models to the ready 23 set (clustering doesn't recognize new addit# ions) have to delete all template sequences and only paste in the extra one we want to do here
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
    
    
