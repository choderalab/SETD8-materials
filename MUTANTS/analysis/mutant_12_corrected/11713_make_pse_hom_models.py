run = 12

cmd.load('files/mutant_homology_models/%d.pdb' % run, 'traj-run%d' % run)
cmd.load_traj('mutant_12_corrected_new_frames_traj_noH.pdb', 'traj-run%d' % run, state=1)
cmd.intra_fit('traj-run%d' % run)
cmd.align('traj-run%d' % run, '1zkk')
cmd.hide('everything', 'traj-run%d' % run)
cmd.show('cartoon', 'traj-run%d' % run)
    
# make custom colors
cmd.set_color('a2', [137, 120, 181])
cmd.set_color('b1', [234, 112, 21])
cmd.set_color('c3', [209, 233, 207])
cmd.set_color('d3', [223, 218, 168])
# color all residues
cmd.color('d3', 'all')
# color selected residues
cmd.color('b1', 'resi 146-162')
cmd.color('a2', 'resi 60-94')
cmd.color('c3', 'resi 1-25')

cmd.save('run12_corrected_new_frames.pse')
