import mdtraj as md
import pdbfixer
from simtk.openmm.app import PDBFile, Modeller, ForceField
import os
import numpy as np

if not os.path.exists('temp_structures/'):
    os.mkdir('temp_structures')
if not os.path.exists('no_solvent/'):
    os.mkdir('no_solvent')
if not os.path.exists('with_solvent/'):
    os.mkdir('with_solvent')

#####
# first we will generate structures based off the 2BQZ, then in section 2 off 4IJ8
# then in section 3 off 1ZKK_apo domain chimera
# Sections separated with lines of hashes
#####


#################
###
# 1ZKK_apo based structures - with all combinations as in the 2BQZ set-up
###
# Copied the model.pdb for 1ZKK_apo from Ensembler set-up of p11709,
# plan is here to superpose the Nflank, SET(1), SET-I, SET(2)-2AA part with 2BQZ, and stuck on all the same ligands,
# then all the sam_modeller
###

# loading here Ensembler model - so 1 chain, already protonated
zkk_apo = md.load('files/1ZKK_apo_Ensembler_model_p11709.pdb')

# remove hydrogens from zkk_apo - so that all structures here are set-up by
# protonation in the final state of all ligands
# thinking is protonation states might be different between RUNs, but they will
# be appropriate to the conformational space around the starting structures
# and that's what we're primarily going for

# also remove residue 0 from 1ZKK - the Ensembler runs had one more AA at the N-terminus
# compared to the 2BQZ chain0 sequence with which we are going here

zkk_apo_noH_selection = zkk_apo.top.select('not element H and not resid 0')
zkk_apo = zkk_apo.atom_slice(zkk_apo_noH_selection)
# we save and re-load to reset the resid's (now that we've cut off the first residue)
# the resid of the new first residue remains 1
# this will be confusing later on unless we reset
zkk_apo.save('temp_structures/1ZKK_apo_SET8_noH.pdb')

zkk_apo_set8 = md.load('temp_structures/1ZKK_apo_SET8_noH.pdb')

# superpose the zkk_apo onto bqz
set8 = md.load('temp_structures/SET8_noH.pdb')

# for the superpose we are NOT going to consider: the Cflank (obviously different
# as 1zkk_apo is a chimera wth apo Cflank) - Cflank same as in 11709 is resid >= 145 - so do here resid < 145
# Nflank - just about maybe 10 residues, at the very N-terminus deviate a bit,
# cut off whole Nflank for convenience - that is first 25 residues - do resid > 24
atom_indices = zkk_apo_set8.top.select('resid > 24 and resid < 145')
ref_atom_indices = set8.top.select('resid > 24 and resid < 145')

zkk_apo_set8 = zkk_apo_set8.superpose(set8, atom_indices=atom_indices, ref_atom_indices=ref_atom_indices)

# separate the P and MeP from SET8_P and SET8_MeP
set8_p = md.load('temp_structures/SET8_P_noH.pdb')
set8_mep = md.load('temp_structures/SET8_MeP_noH.pdb')

p_selection = set8_p.top.select('chainid 1')
p = set8_p.atom_slice(p_selection)
mep_selection = set8_mep.top.select('chainid 1')
mep = set8_mep.atom_slice(mep_selection)

# load in sam and sah
sam = md.load('temp_structures/SAM_noH.pdb')
sah = md.load('temp_structures/SAH_noH.pdb')

# stack all structures
# elements are
# zkk_apo_set8, p, mep, sam, sah
# 9 structures at this stage
# zkk_apo_set8 ready

zkk_apo_set8_p = zkk_apo_set8.stack(p)
zkk_apo_set8_mep = zkk_apo_set8.stack(mep)

zkk_apo_set8_sam = zkk_apo_set8.stack(sam)
zkk_apo_set8_sah = zkk_apo_set8.stack(sah)

zkk_apo_set8_p_sam = zkk_apo_set8_p.stack(sam)
zkk_apo_set8_p_sah = zkk_apo_set8_p.stack(sah)

zkk_apo_set8_mep_sam = zkk_apo_set8_mep.stack(sam)
zkk_apo_set8_mep_sah = zkk_apo_set8_mep.stack(sah)

# save the new structures
zkk_apo_set8.save('temp_structures/1ZKK_apo_SET8_noH.pdb')
zkk_apo_set8_p.save('temp_structures/1ZKK_apo_SET8_P_noH.pdb')
zkk_apo_set8_mep.save('temp_structures/1ZKK_apo_SET8_MeP_noH.pdb')

zkk_apo_set8_sam.save('temp_structures/1ZKK_apo_SET8_SAM_noH.pdb')
zkk_apo_set8_sah.save('temp_structures/1ZKK_apo_SET8_SAH_noH.pdb')

zkk_apo_set8_p_sam.save('temp_structures/1ZKK_apo_SET8_P_SAM_noH.pdb')
zkk_apo_set8_p_sah.save('temp_structures/1ZKK_apo_SET8_P_SAH_noH.pdb')

zkk_apo_set8_mep_sam.save('temp_structures/1ZKK_apo_SET8_MeP_SAM_noH.pdb')
zkk_apo_set8_mep_sah.save('temp_structures/1ZKK_apo_SET8_MeP_SAH_noH.pdb')

# protonate all 9 structures (and SAM/SAH) - through simtk.openmm.app.Modeller
# MDTraj topology is no good for creating Modeller objects - gotta load into PDBFile objects
zkk_apo_set8 = PDBFile('temp_structures/1ZKK_apo_SET8_noH.pdb')
zkk_apo_set8_p = PDBFile('temp_structures/1ZKK_apo_SET8_P_noH.pdb')
zkk_apo_set8_mep = PDBFile('temp_structures/1ZKK_apo_SET8_MeP_noH.pdb')

zkk_apo_set8_sam = PDBFile('temp_structures/1ZKK_apo_SET8_SAM_noH.pdb')
zkk_apo_set8_sah = PDBFile('temp_structures/1ZKK_apo_SET8_SAH_noH.pdb')

zkk_apo_set8_p_sam = PDBFile('temp_structures/1ZKK_apo_SET8_P_SAM_noH.pdb')
zkk_apo_set8_p_sah = PDBFile('temp_structures/1ZKK_apo_SET8_P_SAH_noH.pdb')

zkk_apo_set8_mep_sam = PDBFile('temp_structures/1ZKK_apo_SET8_MeP_SAM_noH.pdb')
zkk_apo_set8_mep_sah = PDBFile('temp_structures/1ZKK_apo_SET8_MeP_SAH_noH.pdb')

# create Modeller objects
# P is peptide with deprotonated (LYN) reactive lysine, PH with protonated (LYS)
# we are making both versions of all P containing species - just so we have them
# decision which ones are viable and will be simulated - later
# 11 final structures here after adding PH variant

zkk_apo_set8_modeller = Modeller(zkk_apo_set8.topology, zkk_apo_set8.positions)
zkk_apo_set8_p_modeller = Modeller(zkk_apo_set8_p.topology, zkk_apo_set8_p.positions)
zkk_apo_set8_ph_modeller = Modeller(zkk_apo_set8_p.topology, zkk_apo_set8_p.positions)
zkk_apo_set8_mep_modeller = Modeller(zkk_apo_set8_mep.topology, zkk_apo_set8_mep.positions)

zkk_apo_set8_sam_modeller = Modeller(zkk_apo_set8_sam.topology, zkk_apo_set8_sam.positions)
zkk_apo_set8_sah_modeller = Modeller(zkk_apo_set8_sah.topology, zkk_apo_set8_sah.positions)

zkk_apo_set8_p_sam_modeller = Modeller(zkk_apo_set8_p_sam.topology, zkk_apo_set8_p_sam.positions)
zkk_apo_set8_ph_sam_modeller = Modeller(zkk_apo_set8_p_sam.topology, zkk_apo_set8_p_sam.positions)
zkk_apo_set8_p_sah_modeller = Modeller(zkk_apo_set8_p_sah.topology, zkk_apo_set8_p_sah.positions)
zkk_apo_set8_ph_sah_modeller = Modeller(zkk_apo_set8_p_sah.topology, zkk_apo_set8_p_sah.positions)

zkk_apo_set8_mep_sam_modeller = Modeller(zkk_apo_set8_mep_sam.topology, zkk_apo_set8_mep_sam.positions)
zkk_apo_set8_mep_sah_modeller = Modeller(zkk_apo_set8_mep_sah.topology, zkk_apo_set8_mep_sah.positions)

# load hydrogen.xml definitions
zkk_apo_set8_modeller.loadHydrogenDefinitions('files/hydrogens.xml')
zkk_apo_set8_p_modeller.loadHydrogenDefinitions('files/hydrogens.xml')
zkk_apo_set8_ph_modeller.loadHydrogenDefinitions('files/hydrogens.xml')
zkk_apo_set8_mep_modeller.loadHydrogenDefinitions('files/hydrogens.xml')

zkk_apo_set8_sam_modeller.loadHydrogenDefinitions('files/hydrogens.xml')
zkk_apo_set8_sah_modeller.loadHydrogenDefinitions('files/hydrogens.xml')

zkk_apo_set8_p_sam_modeller.loadHydrogenDefinitions('files/hydrogens.xml')
zkk_apo_set8_ph_sam_modeller.loadHydrogenDefinitions('files/hydrogens.xml')
zkk_apo_set8_p_sah_modeller.loadHydrogenDefinitions('files/hydrogens.xml')
zkk_apo_set8_ph_sah_modeller.loadHydrogenDefinitions('files/hydrogens.xml')

zkk_apo_set8_mep_sam_modeller.loadHydrogenDefinitions('files/hydrogens.xml')
zkk_apo_set8_mep_sah_modeller.loadHydrogenDefinitions('files/hydrogens.xml')

# add hydrogens - create variants list for anything containing P - need neutral reactive lysine on the peptide (LYN)
# and also make PH (with protonated LYS) versions of those by not specifying variants
zkk_apo_set8_modeller.addHydrogens()
zkk_apo_set8_p_variants = [None] * 171
zkk_apo_set8_p_variants[164] = 'LYN'
zkk_apo_set8_p_modeller.addHydrogens(variants=zkk_apo_set8_p_variants)
zkk_apo_set8_ph_modeller.addHydrogens()
zkk_apo_set8_mep_modeller.addHydrogens()

zkk_apo_set8_sam_modeller.addHydrogens()
zkk_apo_set8_sah_modeller.addHydrogens()

zkk_apo_set8_p_sam_or_sah_variants = [None] * 172
zkk_apo_set8_p_sam_or_sah_variants[164] = 'LYN'
zkk_apo_set8_p_sam_modeller.addHydrogens(variants=zkk_apo_set8_p_sam_or_sah_variants)
zkk_apo_set8_ph_sam_modeller.addHydrogens()
zkk_apo_set8_p_sah_modeller.addHydrogens(variants=zkk_apo_set8_p_sam_or_sah_variants)
zkk_apo_set8_ph_sah_modeller.addHydrogens()

zkk_apo_set8_mep_sam_modeller.addHydrogens()
zkk_apo_set8_mep_sah_modeller.addHydrogens()

# save all in no_solvent - 9 structures
PDBFile.writeFile(zkk_apo_set8_modeller.topology, zkk_apo_set8_modeller.positions, open('no_solvent/1ZKK_apo_SET8.pdb', 'w'))
PDBFile.writeFile(zkk_apo_set8_p_modeller.topology, zkk_apo_set8_p_modeller.positions, open('no_solvent/1ZKK_apo_SET8_P.pdb','w'))
PDBFile.writeFile(zkk_apo_set8_ph_modeller.topology, zkk_apo_set8_ph_modeller.positions, open('no_solvent/1ZKK_apo_SET8_PH.pdb','w'))
PDBFile.writeFile(zkk_apo_set8_mep_modeller.topology, zkk_apo_set8_mep_modeller.positions, open('no_solvent/1ZKK_apo_SET8_MeP.pdb', 'w'))

PDBFile.writeFile(zkk_apo_set8_sam_modeller.topology, zkk_apo_set8_sam_modeller.positions, open('no_solvent/1ZKK_apo_SET8_SAM.pdb','w'))
PDBFile.writeFile(zkk_apo_set8_sah_modeller.topology, zkk_apo_set8_sah_modeller.positions, open('no_solvent/1ZKK_apo_SET8_SAH.pdb', 'w'))

PDBFile.writeFile(zkk_apo_set8_p_sam_modeller.topology, zkk_apo_set8_p_sam_modeller.positions, open('no_solvent/1ZKK_apo_SET8_P_SAM.pdb','w'))
PDBFile.writeFile(zkk_apo_set8_ph_sam_modeller.topology, zkk_apo_set8_ph_sam_modeller.positions, open('no_solvent/1ZKK_apo_SET8_PH_SAM.pdb','w'))
PDBFile.writeFile(zkk_apo_set8_p_sah_modeller.topology, zkk_apo_set8_p_sah_modeller.positions, open('no_solvent/1ZKK_apo_SET8_P_SAH.pdb', 'w'))
PDBFile.writeFile(zkk_apo_set8_ph_sah_modeller.topology, zkk_apo_set8_ph_sah_modeller.positions, open('no_solvent/1ZKK_apo_SET8_PH_SAH.pdb', 'w'))

PDBFile.writeFile(zkk_apo_set8_mep_sam_modeller.topology, zkk_apo_set8_mep_sam_modeller.positions, open('no_solvent/1ZKK_apo_SET8_MeP_SAM.pdb','w'))
PDBFile.writeFile(zkk_apo_set8_mep_sah_modeller.topology, zkk_apo_set8_mep_sah_modeller.positions, open('no_solvent/1ZKK_apo_SET8_MeP_SAH.pdb', 'w'))
