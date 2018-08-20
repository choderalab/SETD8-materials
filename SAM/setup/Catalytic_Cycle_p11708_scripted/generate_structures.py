import mdtraj as md
import pdbfixer
from simtk.openmm.app import PDBFile, Modeller, ForceField
import os

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

###
# 2BQZ based structures - with all combinations of ligands
###

# create mdtraj trajectories
bqz = md.load_pdb('files/pdb2bqz.ent')
ij8 = md.load_pdb('files/pdb4ij8.ent')

# create pdbfixer objects
fixer_SET8 = pdbfixer.PDBFixer('files/pdb2bqz.ent')
fixer_P = pdbfixer.PDBFixer('files/pdb2bqz.ent')
fixer_MeP = pdbfixer.PDBFixer('files/pdb2bqz.ent')

# remove all undesired chains
fixer_SET8.removeChains([1,2,3,4,5,6,7,8,9])
fixer_P.removeChains([0,2,3,4,5,6,7,8,9])
fixer_MeP.removeChains([0,2,3,4,5,6,7,8,9])

# generate missingResidues and add SER to add to N-terminus
# for sequence to match that of positions 232-393 - same as in p11707 and p11709
# (the Ensembler set-ups apo and domain swaps)
fixer_SET8.findMissingResidues()
fixer_SET8.missingResidues = {(0,0): ['SER']}

# do addmissingAtoms which will add the SER at 0 position
fixer_SET8.findMissingAtoms()
fixer_SET8.addMissingAtoms()

# replace MLZ with LYS for SET8_P
fixer_P.findNonstandardResidues()
fixer_P.replaceNonstandardResidues()

# save
PDBFile.writeFile(fixer_SET8.topology, fixer_SET8.positions, open('temp_structures/SET8_noH.pdb', 'w'))
PDBFile.writeFile(fixer_P.topology, fixer_P.positions, open('temp_structures/P_noH.pdb', 'w'))
PDBFile.writeFile(fixer_MeP.topology, fixer_MeP.positions, open('temp_structures/MeP_noH.pdb', 'w'))

# the SAH in 2BQZ is missing OXT atom - fix by taking SAM from 4IJ8 (identical conformations) and superposing - that will be SAM, remove methyl and change name to SAH - that will be SAH
ij8_sam_selection = ij8.top.select('chainid 3 and resname SAM')
bqz_sah_selection = bqz.top.select('chainid 4')
sam = ij8.atom_slice(ij8_sam_selection)
sah = bqz.atom_slice(bqz_sah_selection)

# superpose - remove from selection atom CE (of course not present in SAH) and OXT (missing in the 2BQZ SAH)
sam_superpose_selection = sam.top.select('not name CE and not name OXT')
sah_superpose_selection = sah.top.select('all')
sam = sam.superpose(sah, atom_indices=sam_superpose_selection, ref_atom_indices=sah_superpose_selection)
new_sah_selection = sam.top.select('not name CE')
sah = sam.atom_slice(new_sah_selection)
# change name of the residue in sah from SAM to SAH
sah.top.residue(0).name = 'SAH'

# save
sam.save('temp_structures/SAM_noH.pdb')
sah.save('temp_structures/SAH_noH.pdb')

# MDTraj load remaining components - set8, p, mep - sam and sah already loaded
set8 = md.load('temp_structures/SET8_noH.pdb')
p = md.load('temp_structures/P_noH.pdb')
mep = md.load('temp_structures/MeP_noH.pdb')

# make 8 further structures
set8_p = set8.stack(p)
set8_mep = set8.stack(mep)

set8_sam = set8.stack(sam)
set8_sah = set8.stack(sah)

set8_p_sam = set8_p.stack(sam)
set8_p_sah = set8_p.stack(sah)

set8_mep_sam = set8_mep.stack(sam)
set8_mep_sah = set8_mep.stack(sah)

# save the new structures
set8_p.save('temp_structures/SET8_P_noH.pdb')
set8_mep.save('temp_structures/SET8_MeP_noH.pdb')

set8_sam.save('temp_structures/SET8_SAM_noH.pdb')
set8_sah.save('temp_structures/SET8_SAH_noH.pdb')

set8_p_sam.save('temp_structures/SET8_P_SAM_noH.pdb')
set8_p_sah.save('temp_structures/SET8_P_SAH_noH.pdb')

set8_mep_sam.save('temp_structures/SET8_MeP_SAM_noH.pdb')
set8_mep_sah.save('temp_structures/SET8_MeP_SAH_noH.pdb')

# protonate all 9 final structures (and SAM/SAH - need model for antechamber) - through simtk.openmm.app.Modeller
# MDTraj topology is no good for creating Modeller objects - gotta load into PDBFile objects
set8 = PDBFile('temp_structures/SET8_noH.pdb')
set8_p = PDBFile('temp_structures/SET8_P_noH.pdb')
set8_mep = PDBFile('temp_structures/SET8_MeP_noH.pdb')

set8_sam = PDBFile('temp_structures/SET8_SAM_noH.pdb')
set8_sah = PDBFile('temp_structures/SET8_SAH_noH.pdb')

set8_p_sam = PDBFile('temp_structures/SET8_P_SAM_noH.pdb')
set8_p_sah = PDBFile('temp_structures/SET8_P_SAH_noH.pdb')

set8_mep_sam = PDBFile('temp_structures/SET8_MeP_SAM_noH.pdb')
set8_mep_sah = PDBFile('temp_structures/SET8_MeP_SAH_noH.pdb')

sam = PDBFile('temp_structures/SAM_noH.pdb')
sah = PDBFile('temp_structures/SAH_noH.pdb')

# create Modeller objects
# we're going to protonate SAM and SAH alone too - we need those for antechamber
# P is peptide with deprotonated (LYN) reactive lysine, PH with protonated (LYS)
# we are making both versions of all P containing species - just so we have them
# decision which ones are viable and will be simulated - later

set8_modeller = Modeller(set8.topology, set8.positions)
set8_p_modeller = Modeller(set8_p.topology, set8_p.positions)
set8_ph_modeller = Modeller(set8_p.topology, set8_p.positions)
set8_mep_modeller = Modeller(set8_mep.topology, set8_mep.positions)

set8_sam_modeller = Modeller(set8_sam.topology, set8_sam.positions)
set8_sah_modeller = Modeller(set8_sah.topology, set8_sah.positions)

set8_p_sam_modeller = Modeller(set8_p_sam.topology, set8_p_sam.positions)
set8_ph_sam_modeller = Modeller(set8_p_sam.topology, set8_p_sam.positions)
set8_p_sah_modeller = Modeller(set8_p_sah.topology, set8_p_sah.positions)
set8_ph_sah_modeller = Modeller(set8_p_sah.topology, set8_p_sah.positions)

set8_mep_sam_modeller = Modeller(set8_mep_sam.topology, set8_mep_sam.positions)
set8_mep_sah_modeller = Modeller(set8_mep_sah.topology, set8_mep_sah.positions)

sam_modeller = Modeller(sam.topology, sam.positions)
sah_modeller = Modeller(sah.topology, sah.positions)

# load hydrogen.xml definitions
set8_modeller.loadHydrogenDefinitions('files/hydrogens.xml')
set8_p_modeller.loadHydrogenDefinitions('files/hydrogens.xml')
set8_ph_modeller.loadHydrogenDefinitions('files/hydrogens.xml')
set8_mep_modeller.loadHydrogenDefinitions('files/hydrogens.xml')

set8_sam_modeller.loadHydrogenDefinitions('files/hydrogens.xml')
set8_sah_modeller.loadHydrogenDefinitions('files/hydrogens.xml')

set8_p_sam_modeller.loadHydrogenDefinitions('files/hydrogens.xml')
set8_ph_sam_modeller.loadHydrogenDefinitions('files/hydrogens.xml')
set8_p_sah_modeller.loadHydrogenDefinitions('files/hydrogens.xml')
set8_ph_sah_modeller.loadHydrogenDefinitions('files/hydrogens.xml')

set8_mep_sam_modeller.loadHydrogenDefinitions('files/hydrogens.xml')
set8_mep_sah_modeller.loadHydrogenDefinitions('files/hydrogens.xml')

sam_modeller.loadHydrogenDefinitions('files/hydrogens.xml')
sah_modeller.loadHydrogenDefinitions('files/hydrogens.xml')

# add hydrogens - create variants list for anything containing P - need neutral reactive lysine on the peptide (LYN)
# and also make PH (with protonated LYS) versions of those by not specifying variants
set8_modeller.addHydrogens()
set8_p_variants = [None] * 172
set8_p_variants[165] = 'LYN'
set8_p_modeller.addHydrogens(variants=set8_p_variants)
set8_ph_modeller.addHydrogens()
set8_mep_modeller.addHydrogens()

set8_sam_modeller.addHydrogens()
set8_sah_modeller.addHydrogens()

set8_p_sam_or_sah_variants = [None] * 173
set8_p_sam_or_sah_variants[165] = 'LYN'
set8_p_sam_modeller.addHydrogens(variants=set8_p_sam_or_sah_variants)
set8_ph_sam_modeller.addHydrogens()
set8_p_sah_modeller.addHydrogens(variants=set8_p_sam_or_sah_variants)
set8_ph_sah_modeller.addHydrogens()

set8_mep_sam_modeller.addHydrogens()
set8_mep_sah_modeller.addHydrogens()

sam_modeller.addHydrogens()
sah_modeller.addHydrogens()

# save all in no_solvent - 12 structures
PDBFile.writeFile(set8_modeller.topology, set8_modeller.positions, open('no_solvent/SET8.pdb', 'w'))
PDBFile.writeFile(set8_p_modeller.topology, set8_p_modeller.positions, open('no_solvent/SET8_P.pdb','w'))
PDBFile.writeFile(set8_ph_modeller.topology, set8_ph_modeller.positions, open('no_solvent/SET8_PH.pdb','w'))
PDBFile.writeFile(set8_mep_modeller.topology, set8_mep_modeller.positions, open('no_solvent/SET8_MeP.pdb', 'w'))

PDBFile.writeFile(set8_sam_modeller.topology, set8_sam_modeller.positions, open('no_solvent/SET8_SAM.pdb','w'))
PDBFile.writeFile(set8_sah_modeller.topology, set8_sah_modeller.positions, open('no_solvent/SET8_SAH.pdb', 'w'))

PDBFile.writeFile(set8_p_sam_modeller.topology, set8_p_sam_modeller.positions, open('no_solvent/SET8_P_SAM.pdb','w'))
PDBFile.writeFile(set8_ph_sam_modeller.topology, set8_ph_sam_modeller.positions, open('no_solvent/SET8_PH_SAM.pdb','w'))
PDBFile.writeFile(set8_p_sah_modeller.topology, set8_p_sah_modeller.positions, open('no_solvent/SET8_P_SAH.pdb', 'w'))
PDBFile.writeFile(set8_ph_sah_modeller.topology, set8_ph_sah_modeller.positions, open('no_solvent/SET8_PH_SAH.pdb', 'w'))

PDBFile.writeFile(set8_mep_sam_modeller.topology, set8_mep_sam_modeller.positions, open('no_solvent/SET8_MeP_SAM.pdb','w'))
PDBFile.writeFile(set8_mep_sah_modeller.topology, set8_mep_sah_modeller.positions, open('no_solvent/SET8_MeP_SAH.pdb', 'w'))

PDBFile.writeFile(sam_modeller.topology, sam_modeller.positions, open('no_solvent/SAM.pdb','w'))
PDBFile.writeFile(sah_modeller.topology, sah_modeller.positions, open('no_solvent/SAH.pdb', 'w'))

#################
###
# 4IJ8 based structures - with SAM and SAH
###
# 4ij8 has the right length in chain0, except for one mutation - 0-indexed position 100 - change S to C
# do not add any missing residues - would add 4 unneeded residues at N-terminus
# 
# discovered that the number of waters added for the 4IJ8 systems was bigger than for the 2BQZ and 1ZKK systems
# because the protein is not rotated by openmm.app.Modeller to mimize volume
# hence even though we don't need to superpose the 4IJ8 onto 2BQZ, because 4IJ8 is where SAM originally comes from
# we will to get a similar number of waters throughout
###

fixer_4ij8_protein = pdbfixer.PDBFixer('files/pdb4ij8.ent')
fixer_4ij8_protein.removeChains([1,2,3,4,5,6])
fixer_4ij8_protein.findMissingResidues()
# only add an S to have same sequence as p11707 and p11709
fixer_4ij8_protein.missingResidues = {(0,0): ['SER']}
fixer_4ij8_protein.findNonstandardResidues()
residues = list(fixer_4ij8_protein.topology.residues())
fixer_4ij8_protein.nonstandardResidues = [(residues[110], 'CYS')]
fixer_4ij8_protein.replaceNonstandardResidues()
fixer_4ij8_protein.findMissingAtoms()
fixer_4ij8_protein.addMissingAtoms()
PDBFile.writeFile(fixer_4ij8_protein.topology, fixer_4ij8_protein.positions, open('temp_structures/4IJ8_SET8_noH.pdb', 'w'))

# make sam and sah
ij8 = md.load_pdb('files/pdb4ij8.ent')
ij8_sam_selection = ij8.top.select('chainid 3 and resname SAM')
sam = ij8.atom_slice(ij8_sam_selection)
new_sah_selection = sam.top.select('not name CE')
sah = sam.atom_slice(new_sah_selection)
sah.top.residue(0).name = 'SAH'

# load in ij8_set8 and stack with sam,sah
ij8_set8 = md.load('temp_structures/4IJ8_SET8_noH.pdb')
ij8_set8_sam = ij8_set8.stack(sam)
ij8_set8_sah = ij8_set8.stack(sah)

# superpose onto 2BQZ - to get same rotation and hence similar number of waters
# superpose BY LIGANDS 
set8_sam = md.load('temp_structures/SET8_SAM_noH.pdb')
set8_sah = md.load('temp_structures/SET8_SAH_noH.pdb')

ij8_set8_sam_selection = ij8_set8_sam.top.select('chainid 1')
ij8_set8_sah_selection = ij8_set8_sah.top.select('chainid 1')
set8_sam_selection = set8_sam.top.select('chainid 1')
set8_sah_selection = set8_sah.top.select('chainid 1')

ij8_set8_sam = ij8_set8_sam.superpose(set8_sam, atom_indices=ij8_set8_sam_selection, ref_atom_indices=set8_sam_selection)
ij8_set8_sah = ij8_set8_sah.superpose(set8_sah, atom_indices=ij8_set8_sah_selection, ref_atom_indices=set8_sah_selection)

# now supepose the ij8_set8 onto protein
ij8_set8_selection = ij8_set8.top.select('all')
ij8_set8_sam_nosam_selection = ij8_set8_sam.top.select('chainid 0')
ij8_set8 = ij8_set8.superpose(ij8_set8_sam, atom_indices=ij8_set8_selection, ref_atom_indices=ij8_set8_sam_nosam_selection)

ij8_set8.save('temp_structures/4IJ8_SET8_noH.pdb')
ij8_set8_sam.save('temp_structures/4IJ8_SET8_SAM_noH.pdb')
ij8_set8_sah.save('temp_structures/4IJ8_SET8_SAH_noH.pdb')

ij8_set8 = PDBFile('temp_structures/4IJ8_SET8_noH.pdb')
ij8_set8_sam = PDBFile('temp_structures/4IJ8_SET8_SAM_noH.pdb')
ij8_set8_sah = PDBFile('temp_structures/4IJ8_SET8_SAH_noH.pdb')

ij8_set8_modeller = Modeller(ij8_set8.topology, ij8_set8.positions)
ij8_set8_sam_modeller = Modeller(ij8_set8_sam.topology, ij8_set8_sam.positions)
ij8_set8_sah_modeller = Modeller(ij8_set8_sah.topology, ij8_set8_sah.positions)

ij8_set8_modeller.loadHydrogenDefinitions('files/hydrogens.xml')
ij8_set8_sam_modeller.loadHydrogenDefinitions('files/hydrogens.xml')
ij8_set8_sah_modeller.loadHydrogenDefinitions('files/hydrogens.xml')

ij8_set8_modeller.addHydrogens()
ij8_set8_sam_modeller.addHydrogens()
ij8_set8_sah_modeller.addHydrogens()

PDBFile.writeFile(ij8_set8_modeller.topology, ij8_set8_modeller.positions, open('no_solvent/4IJ8_SET8.pdb','w'))
PDBFile.writeFile(ij8_set8_sam_modeller.topology, ij8_set8_sam_modeller.positions, open('no_solvent/4IJ8_SET8_SAM.pdb','w'))
PDBFile.writeFile(ij8_set8_sah_modeller.topology, ij8_set8_sah_modeller.positions, open('no_solvent/4IJ8_SET8_SAH.pdb', 'w'))

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

zkk_apo_noH_selection = zkk_apo.top.select('not element H')
zkk_apo_set8 = zkk_apo.atom_slice(zkk_apo_noH_selection)

# superpose the zkk_apo onto bqz
set8 = md.load('temp_structures/SET8_noH.pdb')
# for the superpose we are NOT going to consider: the Cflank (obviously different
# as 1zkk_apo is a chimera wth apo Cflank) - Cflank same as in 11709 is resid >= 145 - so do here resid < 145
# Nflank - just about maybe 10 residues, at the very N-terminus deviate a bit,
# cut off whole Nflank for convenience - that is first 25 residues - do resid > 24
atom_indices = zkk_apo_set8.top.select('resid > 24 and resid < 145')
ref_atom_indices = set8.top.select('resid > 24 and resid < 145')

zkk_apo_set8 = zkk_apo_set8.superpose(set8, atom_indices=atom_indices, ref_atom_indices=ref_atom_indices)

# load in p, mep,  sam and sah
p = md.load('temp_structures/P_noH.pdb')
mep = md.load('temp_structures/MeP_noH.pdb')
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
# 12 final structures here after adding PH variant

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
zkk_apo_set8_p_variants = [None] * 172
zkk_apo_set8_p_variants[165] = 'LYN'
zkk_apo_set8_p_modeller.addHydrogens(variants=zkk_apo_set8_p_variants)
zkk_apo_set8_ph_modeller.addHydrogens()
zkk_apo_set8_mep_modeller.addHydrogens()

zkk_apo_set8_sam_modeller.addHydrogens()
zkk_apo_set8_sah_modeller.addHydrogens()

zkk_apo_set8_p_sam_or_sah_variants = [None] * 173
zkk_apo_set8_p_sam_or_sah_variants[165] = 'LYN'
zkk_apo_set8_p_sam_modeller.addHydrogens(variants=zkk_apo_set8_p_sam_or_sah_variants)
zkk_apo_set8_ph_sam_modeller.addHydrogens()
zkk_apo_set8_p_sah_modeller.addHydrogens(variants=zkk_apo_set8_p_sam_or_sah_variants)
zkk_apo_set8_ph_sah_modeller.addHydrogens()

zkk_apo_set8_mep_sam_modeller.addHydrogens()
zkk_apo_set8_mep_sah_modeller.addHydrogens()

# save all in no_solvent - 12 structures
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
