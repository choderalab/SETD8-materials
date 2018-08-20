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

# create mdtraj trajectories
bqz = md.load_pdb('files/pdb2bqz.ent')
ij8 = md.load_pdb('files/pdb4ij8.ent')

# create pdbfixer objects
fixer_SET8 = pdbfixer.PDBFixer('files/pdb2bqz.ent')
fixer_SET8_P = pdbfixer.PDBFixer('files/pdb2bqz.ent')
fixer_SET8_MeP = pdbfixer.PDBFixer('files/pdb2bqz.ent')

# remove all undesired chains
fixer_SET8.removeChains([1,2,3,4,5,6,7,8,9])
fixer_SET8_P.removeChains([2,3,4,5,6,7,8,9])
fixer_SET8_MeP.removeChains([2,3,4,5,6,7,8,9])

# replace MLZ with LYS for SET8_P
fixer_SET8_P.findNonstandardResidues()
fixer_SET8_P.replaceNonstandardResidues()

# save
PDBFile.writeFile(fixer_SET8.topology, fixer_SET8.positions, open('temp_structures/SET8_noH.pdb', 'w'))
PDBFile.writeFile(fixer_SET8_P.topology, fixer_SET8_P.positions, open('temp_structures/SET8_P_noH.pdb', 'w'))
PDBFile.writeFile(fixer_SET8_MeP.topology, fixer_SET8_MeP.positions, open('temp_structures/SET8_MeP_noH.pdb', 'w'))

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

# MDTraj load remaining components - set8, set8_p, set8_mep - sam and sah already loaded
set8 = md.load('temp_structures/SET8_noH.pdb')
set8_p = md.load('temp_structures/SET8_P_noH.pdb')
set8_mep = md.load('temp_structures/SET8_MeP_noH.pdb')

# make 6 further structures
set8_sam = set8.stack(sam)
set8_sah = set8.stack(sah)

set8_p_sam = set8_p.stack(sam)
set8_p_sah = set8_p.stack(sah)

set8_mep_sam = set8_mep.stack(sam)
set8_mep_sah = set8_mep.stack(sah)

# save the new structures
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
set8_modeller = Modeller(set8.topology, set8.positions)
set8_p_modeller = Modeller(set8_p.topology, set8_p.positions)
set8_mep_modeller = Modeller(set8_mep.topology, set8_mep.positions)

set8_sam_modeller = Modeller(set8_sam.topology, set8_sam.positions)
set8_sah_modeller = Modeller(set8_sah.topology, set8_sah.positions)

set8_p_sam_modeller = Modeller(set8_p_sam.topology, set8_p_sam.positions)
set8_p_sah_modeller = Modeller(set8_p_sah.topology, set8_p_sah.positions)

set8_mep_sam_modeller = Modeller(set8_mep_sam.topology, set8_mep_sam.positions)
set8_mep_sah_modeller = Modeller(set8_mep_sah.topology, set8_mep_sah.positions)

sam_modeller = Modeller(sam.topology, sam.positions)
sah_modeller = Modeller(sah.topology, sah.positions)

# load hydrogen.xml definitions
set8_modeller.loadHydrogenDefinitions('files/hydrogens.xml')
set8_p_modeller.loadHydrogenDefinitions('files/hydrogens.xml') 
set8_mep_modeller.loadHydrogenDefinitions('files/hydrogens.xml')

set8_sam_modeller.loadHydrogenDefinitions('files/hydrogens.xml')
set8_sah_modeller.loadHydrogenDefinitions('files/hydrogens.xml')

set8_p_sam_modeller.loadHydrogenDefinitions('files/hydrogens.xml')
set8_p_sah_modeller.loadHydrogenDefinitions('files/hydrogens.xml')

set8_mep_sam_modeller.loadHydrogenDefinitions('files/hydrogens.xml')
set8_mep_sah_modeller.loadHydrogenDefinitions('files/hydrogens.xml')

sam_modeller.loadHydrogenDefinitions('files/hydrogens.xml')
sah_modeller.loadHydrogenDefinitions('files/hydrogens.xml')

# add hydrogens - create variants list for anything containing P - need neutral reactive lysine on the peptide (LYN)
set8_modeller.addHydrogens()
set8_p_variants = [None] * 171
set8_p_variants[164] = 'LYN'
set8_p_modeller.addHydrogens(variants=set8_p_variants)
set8_mep_modeller.addHydrogens()

set8_sam_modeller.addHydrogens()
set8_sah_modeller.addHydrogens()

set8_p_sam_or_sah_variants = [None] * 172
set8_p_sam_or_sah_variants[164] = 'LYN'
set8_p_sam_modeller.addHydrogens(variants=set8_p_sam_or_sah_variants)
set8_p_sah_modeller.addHydrogens(variants=set8_p_sam_or_sah_variants)

set8_mep_sam_modeller.addHydrogens()
set8_mep_sah_modeller.addHydrogens()

sam_modeller.addHydrogens()
sah_modeller.addHydrogens()

# save all in no_solvent - 9 structures
PDBFile.writeFile(set8_modeller.topology, set8_modeller.positions, open('no_solvent/SET8.pdb', 'w'))
PDBFile.writeFile(set8_p_modeller.topology, set8_p_modeller.positions, open('no_solvent/SET8_P.pdb','w'))
PDBFile.writeFile(set8_mep_modeller.topology, set8_mep_modeller.positions, open('no_solvent/SET8_MeP.pdb', 'w'))

PDBFile.writeFile(set8_sam_modeller.topology, set8_sam_modeller.positions, open('no_solvent/SET8_SAM.pdb','w'))
PDBFile.writeFile(set8_sah_modeller.topology, set8_sah_modeller.positions, open('no_solvent/SET8_SAH.pdb', 'w'))

PDBFile.writeFile(set8_p_sam_modeller.topology, set8_p_sam_modeller.positions, open('no_solvent/SET8_P_SAM.pdb','w'))
PDBFile.writeFile(set8_p_sah_modeller.topology, set8_p_sah_modeller.positions, open('no_solvent/SET8_P_SAH.pdb', 'w'))

PDBFile.writeFile(set8_mep_sam_modeller.topology, set8_mep_sam_modeller.positions, open('no_solvent/SET8_MeP_SAM.pdb','w'))
PDBFile.writeFile(set8_mep_sah_modeller.topology, set8_mep_sah_modeller.positions, open('no_solvent/SET8_MeP_SAH.pdb', 'w'))

PDBFile.writeFile(sam_modeller.topology, sam_modeller.positions, open('no_solvent/SAM.pdb','w'))
PDBFile.writeFile(sah_modeller.topology, sah_modeller.positions, open('no_solvent/SAH.pdb', 'w'))
