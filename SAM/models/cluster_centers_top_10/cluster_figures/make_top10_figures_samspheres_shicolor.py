#cmd.load('/Users/rafalpwiewiora/repos/MSM/sasa/models/1zkk_a_oriented.pdb', '1zkk')
#cmd.load('/Users/rafalpwiewiora/repos/MSM/sasa/models/Apo_a_oriented.pdb', 'Apo')
#cmd.align('Apo', '1zkk')

#cmd.hide('everything', 'all')
#cmd.show('cartoon', 'all')
#cmd.color('white', 'all')

#cmd.color('black', 'Apo and (resi 60-94 or resi 145-162)')
#cmd.color('red', '1zkk and (resi 60-94 or resi 145-162)')

import glob

for file_ in glob.glob('files/top10/*.dcd'):
    i = int(file_.split('/')[-1][:-4])
    cmd.load('files/4IJ8_SET8_SAM.pdb', '%d' % i)
    cmd.load_traj('files/top10/%d.dcd' % i, '%d' % i, state=1)
    cmd.intra_fit('%d' % i)
    cmd.align('%d' % i, '1zkk')
    cmd.hide('everything', '%d' % i)
    cmd.show('cartoon', '%d' % i)
    cmd.show('spheres', '%d and chain B' % i)
    # make custom colors
    cmd.set_color('a2', [137, 120, 181])
    cmd.set_color('b1', [234, 112, 21])
    cmd.set_color('c3', [209, 233, 207])
    cmd.set_color('d3', [223, 218, 168])
    # color all residues
    cmd.color('d3', '%d' % i)
    # color selected residues
    cmd.color('b1', '%d and resi 146-162' % i)
    cmd.color('a2', '%d and resi 60-94' % i)
    cmd.color('c3', '%d and resi 1-25' % i)
    cmd.color('cyan', '%d and chain B' % i)
    cmd.split_states('%d' % i)
    cmd.delete('%d' % i)
    cmd.set('cartoon_transparency', 0.8)
    cmd.set('sphere_transparency', 0.8)
    cmd.set('cartoon_transparency', 0, '1zkk')
    cmd.set('cartoon_transparency', 0, 'apo')
    cmd.set('cartoon_transparency', 0, '%d_0001' % i)
    cmd.set('sphere_transparency', 0, '%d_0001' % i)
    #cmd.zoom()
    cmd.draw(2400)
    cmd.png('top10_figures/%d.png' % i)
    cmd.delete('%d_0001' % i)
    cmd.delete('%d_0002' % i)
    cmd.delete('%d_0003' % i)
    cmd.delete('%d_0004' % i)
    cmd.delete('%d_0005' % i)
    cmd.delete('%d_0006' % i)
    cmd.delete('%d_0007' % i)
    cmd.delete('%d_0008' % i)
    cmd.delete('%d_0009' % i)
    cmd.delete('%d_0010' % i)
