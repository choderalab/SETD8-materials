import glob

traj = glob.glob('*.pdb')[0]

cmd.load(traj)
cmd.intra_fit(traj[:-4])
x = glob.glob('*.pdb')[0][:-4]

for traj in glob.glob('*.pdb')[1:]:
    cmd.load(traj)
    cmd.intra_fit(traj[:-4])
    cmd.align(traj[:-4], x)

cmd.hide('everything', 'all')
cmd.show('cartoon', 'all')

cmd.zoom()
cmd.dss()

cmd.set_color('a2', [137, 120, 181])
cmd.set_color('b1', [234, 112, 21])
cmd.set_color('c3', [209, 233, 207])
cmd.set_color('d3', [223, 218, 168])

cmd.color('d3')
cmd.color('b1', 'resi 146-162')
cmd.color('a2', 'resi 60-94')
cmd.color('c3', 'resi 1-25')
