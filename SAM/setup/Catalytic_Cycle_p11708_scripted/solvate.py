from simtk.openmm.app import PDBFile, Modeller, ForceField
from simtk import unit as u

ff_sam = ForceField('amber99sbildn.xml', 'tip3p.xml', 'parameters/gaff.xml', 'parameters/ffptm.xml', 'parameters/SAM.xml')
ff_sah = ForceField('amber99sbildn.xml', 'tip3p.xml', 'parameters/gaff.xml', 'parameters/ffptm.xml', 'parameters/SAH.xml')

#####
# first we will generate structures based off the 2BQZ, then in section 2 off 4IJ8
# then in section 3 off 1ZKK_apo domain chimera
# Sections separated with lines of hashes
#####

###
# 2BQZ based structures - with all combinations of ligands
###

# load all 11 structurres in
set8 = PDBFile('no_solvent/SET8.pdb')
set8_p = PDBFile('no_solvent/SET8_P.pdb')
set8_ph = PDBFile('no_solvent/SET8_PH.pdb')
set8_mep = PDBFile('no_solvent/SET8_MeP.pdb')

set8_sam = PDBFile('no_solvent/SET8_SAM.pdb')
set8_sah = PDBFile('no_solvent/SET8_SAH.pdb')

set8_p_sam = PDBFile('no_solvent/SET8_P_SAM.pdb')
set8_ph_sam = PDBFile('no_solvent/SET8_PH_SAM.pdb')
set8_p_sah = PDBFile('no_solvent/SET8_P_SAH.pdb')
set8_ph_sah = PDBFile('no_solvent/SET8_PH_SAH.pdb')

set8_mep_sam = PDBFile('no_solvent/SET8_MeP_SAM.pdb')
set8_mep_sah = PDBFile('no_solvent/SET8_MeP_SAH.pdb')

# create Modeller objects
set8_modeller = Modeller(set8.topology, set8.positions)
set8_p_modeller = Modeller(set8_p.topology, set8_p.positions)
set8_ph_modeller = Modeller(set8_ph.topology, set8_ph.positions)
set8_mep_modeller = Modeller(set8_mep.topology, set8_mep.positions)

set8_sam_modeller = Modeller(set8_sam.topology, set8_sam.positions)
set8_sah_modeller = Modeller(set8_sah.topology, set8_sah.positions)

set8_p_sam_modeller = Modeller(set8_p_sam.topology, set8_p_sam.positions)
set8_ph_sam_modeller = Modeller(set8_ph_sam.topology, set8_ph_sam.positions)
set8_p_sah_modeller = Modeller(set8_p_sah.topology, set8_p_sah.positions)
set8_ph_sah_modeller = Modeller(set8_ph_sah.topology, set8_ph_sah.positions)

set8_mep_sam_modeller = Modeller(set8_mep_sam.topology, set8_mep_sam.positions)
set8_mep_sah_modeller = Modeller(set8_mep_sah.topology, set8_mep_sah.positions)

# solvate
set8_modeller.addSolvent(ff_sam, padding=1*u.nanometers)
set8_p_modeller.addSolvent(ff_sam, padding=1*u.nanometers)
set8_ph_modeller.addSolvent(ff_sam, padding=1*u.nanometers)
set8_mep_modeller.addSolvent(ff_sam, padding=1*u.nanometers)

set8_sam_modeller.addSolvent(ff_sam, padding=1*u.nanometers)
set8_sah_modeller.addSolvent(ff_sah, padding=1*u.nanometers)

set8_p_sam_modeller.addSolvent(ff_sam, padding=1*u.nanometers)
set8_ph_sam_modeller.addSolvent(ff_sam, padding=1*u.nanometers)
set8_p_sah_modeller.addSolvent(ff_sah, padding=1*u.nanometers)
set8_ph_sah_modeller.addSolvent(ff_sah, padding=1*u.nanometers)

set8_mep_sam_modeller.addSolvent(ff_sam, padding=1*u.nanometers)
set8_mep_sah_modeller.addSolvent(ff_sah, padding=1*u.nanometers)

# save all in with_solvent
PDBFile.writeFile(set8_modeller.topology, set8_modeller.positions, open('with_solvent/SET8.pdb', 'w'))
PDBFile.writeFile(set8_p_modeller.topology, set8_p_modeller.positions, open('with_solvent/SET8_P.pdb','w'))
PDBFile.writeFile(set8_ph_modeller.topology, set8_ph_modeller.positions, open('with_solvent/SET8_PH.pdb','w'))
PDBFile.writeFile(set8_mep_modeller.topology, set8_mep_modeller.positions, open('with_solvent/SET8_MeP.pdb', 'w'))

PDBFile.writeFile(set8_sam_modeller.topology, set8_sam_modeller.positions, open('with_solvent/SET8_SAM.pdb','w'))
PDBFile.writeFile(set8_sah_modeller.topology, set8_sah_modeller.positions, open('with_solvent/SET8_SAH.pdb', 'w'))

PDBFile.writeFile(set8_p_sam_modeller.topology, set8_p_sam_modeller.positions, open('with_solvent/SET8_P_SAM.pdb','w'))
PDBFile.writeFile(set8_ph_sam_modeller.topology, set8_ph_sam_modeller.positions, open('with_solvent/SET8_PH_SAM.pdb','w'))
PDBFile.writeFile(set8_p_sah_modeller.topology, set8_p_sah_modeller.positions, open('with_solvent/SET8_P_SAH.pdb', 'w'))
PDBFile.writeFile(set8_ph_sah_modeller.topology, set8_ph_sah_modeller.positions, open('with_solvent/SET8_PH_SAH.pdb', 'w'))

PDBFile.writeFile(set8_mep_sam_modeller.topology, set8_mep_sam_modeller.positions, open('with_solvent/SET8_MeP_SAM.pdb','w'))
PDBFile.writeFile(set8_mep_sah_modeller.topology, set8_mep_sah_modeller.positions, open('with_solvent/SET8_MeP_SAH.pdb', 'w'))

#################
###
# 4IJ8 based structures - with SAM and SAH
###

ij8_set8 = PDBFile('no_solvent/4IJ8_SET8.pdb')
ij8_set8_sam = PDBFile('no_solvent/4IJ8_SET8_SAM.pdb')
ij8_set8_sah = PDBFile('no_solvent/4IJ8_SET8_SAH.pdb')

ij8_set8_modeller = Modeller(ij8_set8.topology, ij8_set8.positions)
ij8_set8_sam_modeller = Modeller(ij8_set8_sam.topology, ij8_set8_sam.positions)
ij8_set8_sah_modeller = Modeller(ij8_set8_sah.topology, ij8_set8_sah.positions)

ij8_set8_modeller.addSolvent(ff_sam, padding=1*u.nanometers)
ij8_set8_sam_modeller.addSolvent(ff_sam, padding=1*u.nanometers)
ij8_set8_sah_modeller.addSolvent(ff_sah, padding=1*u.nanometers)

PDBFile.writeFile(ij8_set8_modeller.topology, ij8_set8_modeller.positions, open('with_solvent/4IJ8_SET8.pdb','w'))
PDBFile.writeFile(ij8_set8_sam_modeller.topology, ij8_set8_sam_modeller.positions, open('with_solvent/4IJ8_SET8_SAM.pdb','w'))
PDBFile.writeFile(ij8_set8_sah_modeller.topology, ij8_set8_sah_modeller.positions, open('with_solvent/4IJ8_SET8_SAH.pdb', 'w'))

#################
###
# 1ZKK_apo based structures - with all combinations as in the 2BQZ set-up
###

# load all 11 structures in
zkk_apo_set8 = PDBFile('no_solvent/1ZKK_apo_SET8.pdb')
zkk_apo_set8_p = PDBFile('no_solvent/1ZKK_apo_SET8_P.pdb')
zkk_apo_set8_ph = PDBFile('no_solvent/1ZKK_apo_SET8_PH.pdb')
zkk_apo_set8_mep = PDBFile('no_solvent/1ZKK_apo_SET8_MeP.pdb')

zkk_apo_set8_sam = PDBFile('no_solvent/1ZKK_apo_SET8_SAM.pdb')
zkk_apo_set8_sah = PDBFile('no_solvent/1ZKK_apo_SET8_SAH.pdb')

zkk_apo_set8_p_sam = PDBFile('no_solvent/1ZKK_apo_SET8_P_SAM.pdb')
zkk_apo_set8_ph_sam = PDBFile('no_solvent/1ZKK_apo_SET8_PH_SAM.pdb')
zkk_apo_set8_p_sah = PDBFile('no_solvent/1ZKK_apo_SET8_P_SAH.pdb')
zkk_apo_set8_ph_sah = PDBFile('no_solvent/1ZKK_apo_SET8_PH_SAH.pdb')

zkk_apo_set8_mep_sam = PDBFile('no_solvent/1ZKK_apo_SET8_MeP_SAM.pdb')
zkk_apo_set8_mep_sah = PDBFile('no_solvent/1ZKK_apo_SET8_MeP_SAH.pdb')

# create Modeller objects
zkk_apo_set8_modeller = Modeller(zkk_apo_set8.topology, zkk_apo_set8.positions)
zkk_apo_set8_p_modeller = Modeller(zkk_apo_set8_p.topology, zkk_apo_set8_p.positions)
zkk_apo_set8_ph_modeller = Modeller(zkk_apo_set8_ph.topology, zkk_apo_set8_ph.positions)
zkk_apo_set8_mep_modeller = Modeller(zkk_apo_set8_mep.topology, zkk_apo_set8_mep.positions)

zkk_apo_set8_sam_modeller = Modeller(zkk_apo_set8_sam.topology, zkk_apo_set8_sam.positions)
zkk_apo_set8_sah_modeller = Modeller(zkk_apo_set8_sah.topology, zkk_apo_set8_sah.positions)

zkk_apo_set8_p_sam_modeller = Modeller(zkk_apo_set8_p_sam.topology, zkk_apo_set8_p_sam.positions)
zkk_apo_set8_ph_sam_modeller = Modeller(zkk_apo_set8_ph_sam.topology, zkk_apo_set8_ph_sam.positions)
zkk_apo_set8_p_sah_modeller = Modeller(zkk_apo_set8_p_sah.topology, zkk_apo_set8_p_sah.positions)
zkk_apo_set8_ph_sah_modeller = Modeller(zkk_apo_set8_ph_sah.topology, zkk_apo_set8_ph_sah.positions)

zkk_apo_set8_mep_sam_modeller = Modeller(zkk_apo_set8_mep_sam.topology, zkk_apo_set8_mep_sam.positions)
zkk_apo_set8_mep_sah_modeller = Modeller(zkk_apo_set8_mep_sah.topology, zkk_apo_set8_mep_sah.positions)

# solvate
zkk_apo_set8_modeller.addSolvent(ff_sam, padding=1*u.nanometers)
zkk_apo_set8_p_modeller.addSolvent(ff_sam, padding=1*u.nanometers)
zkk_apo_set8_ph_modeller.addSolvent(ff_sam, padding=1*u.nanometers)
zkk_apo_set8_mep_modeller.addSolvent(ff_sam, padding=1*u.nanometers)

zkk_apo_set8_sam_modeller.addSolvent(ff_sam, padding=1*u.nanometers)
zkk_apo_set8_sah_modeller.addSolvent(ff_sah, padding=1*u.nanometers)

zkk_apo_set8_p_sam_modeller.addSolvent(ff_sam, padding=1*u.nanometers)
zkk_apo_set8_ph_sam_modeller.addSolvent(ff_sam, padding=1*u.nanometers)
zkk_apo_set8_p_sah_modeller.addSolvent(ff_sah, padding=1*u.nanometers)
zkk_apo_set8_ph_sah_modeller.addSolvent(ff_sah, padding=1*u.nanometers)

zkk_apo_set8_mep_sam_modeller.addSolvent(ff_sam, padding=1*u.nanometers)
zkk_apo_set8_mep_sah_modeller.addSolvent(ff_sah, padding=1*u.nanometers)

# save all in with_solvent
PDBFile.writeFile(zkk_apo_set8_modeller.topology, zkk_apo_set8_modeller.positions, open('with_solvent/1ZKK_apo_SET8.pdb', 'w'))
PDBFile.writeFile(zkk_apo_set8_p_modeller.topology, zkk_apo_set8_p_modeller.positions, open('with_solvent/1ZKK_apo_SET8_P.pdb','w'))
PDBFile.writeFile(zkk_apo_set8_ph_modeller.topology, zkk_apo_set8_ph_modeller.positions, open('with_solvent/1ZKK_apo_SET8_PH.pdb','w'))
PDBFile.writeFile(zkk_apo_set8_mep_modeller.topology, zkk_apo_set8_mep_modeller.positions, open('with_solvent/1ZKK_apo_SET8_MeP.pdb', 'w'))

PDBFile.writeFile(zkk_apo_set8_sam_modeller.topology, zkk_apo_set8_sam_modeller.positions, open('with_solvent/1ZKK_apo_SET8_SAM.pdb','w'))
PDBFile.writeFile(zkk_apo_set8_sah_modeller.topology, zkk_apo_set8_sah_modeller.positions, open('with_solvent/1ZKK_apo_SET8_SAH.pdb', 'w'))

PDBFile.writeFile(zkk_apo_set8_p_sam_modeller.topology, zkk_apo_set8_p_sam_modeller.positions, open('with_solvent/1ZKK_apo_SET8_P_SAM.pdb','w'))
PDBFile.writeFile(zkk_apo_set8_ph_sam_modeller.topology, zkk_apo_set8_ph_sam_modeller.positions, open('with_solvent/1ZKK_apo_SET8_PH_SAM.pdb','w'))
PDBFile.writeFile(zkk_apo_set8_p_sah_modeller.topology, zkk_apo_set8_p_sah_modeller.positions, open('with_solvent/1ZKK_apo_SET8_P_SAH.pdb', 'w'))
PDBFile.writeFile(zkk_apo_set8_ph_sah_modeller.topology, zkk_apo_set8_ph_sah_modeller.positions, open('with_solvent/1ZKK_apo_SET8_PH_SAH.pdb', 'w'))

PDBFile.writeFile(zkk_apo_set8_mep_sam_modeller.topology, zkk_apo_set8_mep_sam_modeller.positions, open('with_solvent/1ZKK_apo_SET8_MeP_SAM.pdb','w'))
PDBFile.writeFile(zkk_apo_set8_mep_sah_modeller.topology, zkk_apo_set8_mep_sah_modeller.positions, open('with_solvent/1ZKK_apo_SET8_MeP_SAH.pdb', 'w'))
