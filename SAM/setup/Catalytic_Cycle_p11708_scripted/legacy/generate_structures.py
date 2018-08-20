import mdtraj as md
import pdbfixer
from simtk.openmm.app import PDBFile

bqz = md.load_pdb('files/pdb2bqz.ent')
ij8 = md.load_pdb('files/pdb4ij8.ent')

fixer_SET8 = pdbfixer.PDBFixer('files/pdb2bqz.ent')
fixer_SET8_P = pdbfixer.PDBFixer('files/pdb2bqz.ent')
fixer_SET8_MeP = pdbfixer.PDBFixer('files/pdb2bqz.ent')

fixer_SET8.removeChains([1,2,3,4,5,6,7,8,9])
fixer_SET8_P.removeChains([2,3,4,5,6,7,8,9])
fixer_SET8_MeP.removeChains([2,3,4,5,6,7,8,9])

fixer_SET8_P.findNonstandardResidues()
fixer_SET8_P.replaceNonstandardResidues()

PDBFile.writeFile(fixer_SET8.topology, fixer_SET8.positions, open('structures/SET8_noH.pdb', 'w'))
PDBFile.writeFile(fixer_SET8_P.topology, fixer_SET8_P.positions, open('structures/SET8_P_noH.pdb', 'w'))
PDBFile.writeFile(fixer_SET8_MeP.topology, fixer_SET8_MeP.positions, open('structures/SET8_MeP_noH.pdb', 'w'))

ij8_sam_selection = ij8.top.select('chainid 3 and resname SAM')
bqz_sah_selection = bqz.top.select('chainid 4')
sam = ij8.atom_slice(ij8_sam_selection)
sah = bqz.atom_slice(bqz_sah_selection)

sam_superpose_selection = sam.top.select('not name CE and not name OXT')
sah_superpose_selection = sah.top.select('all')
sam = sam.superpose(sah, atom_indices=sam_superpose_selection, ref_atom_indices=sah_superpose_selection)
new_sah_selection = sam.top.select('not name CE')
sah = sam.atom_slice(new_sah_selection)

sam.save('structures/SAM_noH.pdb')
sah.save('structures/SAH_noH.pdb')
# reduce will not accept files written by mdtraj hence we need to put them through PDBFile
sam = PDBFile('structures/SAM_noH.pdb')
sah = PDBFile('structures/SAH_noH.pdb')
PDBFile.writeFile(sam.topology, sam.positions, open('structures/SAM_noH.pdb', 'w'))
PDBFile.writeFile(sah.topology, sah.positions, open('structures/SAH_noH.pdb', 'w'))
