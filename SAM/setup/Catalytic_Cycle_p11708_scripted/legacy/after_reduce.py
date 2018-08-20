from simtk.openmm.app import PDBFile
import mdtraj as md
import os

# going to add bonds to the SAM and SAH, taking bonds from PDB of those downloaded from Ligand Expo, by comparing atom names (as order is different)
# TODO: can just change to doing this by mdtraj, we are after reduce so we're ok to just be saving out of mdtraj and we'll do that to all final structures for common format

ideal_sam = PDBFile('files/SAM_ideal.pdb')
ideal_sah = PDBFile('files/SAH_ideal.pdb')

sam = PDBFile('structures/SAM.pdb')
sah = PDBFile('structures/SAH.pdb')

bonds_sam = ideal_sam.topology.bonds()
bonds_sah = ideal_sah.topology.bonds()

for bond in bonds_sam:
    x = bond[0]
    y = bond[1]
    for atom in sam.topology.atoms():
        if atom.name == x.name: 
            x_ = atom
        elif atom.name == y.name:
            y_ = atom
    sam.topology.addBond(x_, y_)
    
for bond in bonds_sah:
    x = bond[0]
    y = bond[1]
    for atom in sah.topology.atoms():
        if atom.name == x.name:
            x_ = atom
        elif atom.name == y.name:
            y_ = atom
    sah.topology.addBond(x_, y_)
    
PDBFile.writeFile(sam.topology, sam.positions, open('structures/SAM.pdb', 'w'))
PDBFile.writeFile(sah.topology, sah.positions, open('structures/SAH.pdb', 'w'))

# also deprotonate the reactive lysine in the P (peptide) - remove atom named HZ3 - atom index 2653
set8_p = md.load('structures/SET8_P.pdb')
set8_p = set8_p.atom_slice(set8_p.top.select('not (chainid 1 and resname LYS and name HZ1)'))
set8_p.save('structures/SET8_P.pdb')

# fix up the protonation of MLZ in SET8_MeP
set8_mep = md.load('structures/SET8_MeP.pdb')
set8_mep = set8_mep.atom_slice(set8_mep.top.select('resname MLZ and (name CM or name HCM1 or name HCM2 or name HCM3)'))
set8_mep = set8_p.stack(set8_mep)
set8_mep.save('structures/SET8_MeP.pdb')

# now combine all the needed elements
if not os.path.exists('no_solvent/'):
    os.mkdir('no_solvent')
    
# setd8_p and setd8_mep already loaded
sam = md.load('structures/SAM.pdb')
sah = md.load('structures/SAH.pdb')
set8 = md.load('structures/SET8.pdb')

# generate mizes of SET8, SET8_P, SET8_MeP with SAM and SAH
set8_sam = set8.stack(sam)
set8_sah = set8.stack(sah)

set8_p_sam = set8_p.stack(sam)
set8_p_sah = set8_p.stack(sah)

set8_mep_sam = set8_mep.stack(sam)
set8_mep_sah = set8_mep.stack(sah)

# save all in no_solvent - 9 structures
set8.save('no_solvent/SET8.pdb')
set8_p.save('no_solvent/SET8_P.pdb')
set8_mep.save('no_solvent/SET8_MeP.pdb')

set8_sam.save('no_solvent/SET8_SAM.pdb')
set8_sah.save('no_solvent/SET8_SAH.pdb')

set8_p_sam.save('no_solvent/SET8_P_SAM.pdb')
set8_p_sah.save('no_solvent/SET8_P_SAH.pdb')

set8_mep_sam.save('no_solvent/SET8_MeP_SAM.pdb')
set8_mep_sah.save('no_solvent/SET8_MeP_SAH.pdb')
