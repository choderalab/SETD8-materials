for i in range(10): cmd.load('/Users/rafalpwiewiora/repos/MSM/models/1zkk_a.pdb', '%d' % i)

for i in range(10): cmd.load_traj('%d.pdb' % i, state=1)

for i in range(10): cmd.intra_fit('%d' % i)

for i in range(9): cmd.align('%d' % i, '9')

cmd.set_color('a2', [137, 120, 181])
cmd.set_color('b1', [234, 112, 21])
cmd.set_color('c3', [209, 233, 207])
cmd.set_color('d3', [223, 218, 168])

cmd.hide('everything', 'all')
cmd.show('cartoon', 'all')

cmd.color('d3', 'all')
# color selected residues
cmd.color('b1', 'resi 146-162')
cmd.color('a2', 'resi 60-94')
cmd.color('c3', 'resi 1-25')
cmd.zoom()
