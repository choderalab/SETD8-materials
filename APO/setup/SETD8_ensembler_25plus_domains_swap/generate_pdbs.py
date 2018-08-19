import mdtraj as md

# we are supplying the models produced in project 11707 with Ensembler from respective crystal structures 
# all are of sequence:  SRKSKAELQSEERKRIDELIESGKEEGMKIDLIDGKGRGVIATKQFSRGDFVVEYHGDLIEITDAKKREALYAQDPSTGCYMYYFQYLSKTYCVDATRETNRLGRLINHSKCGNCQTKLHDIDGVPHLILIASRDIAAGEELLYDYGDRSKASIEAHPWLKH
zkk = md.load_pdb('files/SETD8_HUMAN_1ZKK_A/model.pdb')
ij = md.load_pdb('files/SETD8_HUMAN_4IJ8_A/model.pdb')
apo = md.load_pdb('files/SETD8_HUMAN_APO_A/model.pdb')
inhibitor = md.load_pdb('files/SETD8_HUMAN_INHIBITOR_A/model.pdb')
inhibitor2 = md.load_pdb('files/SETD8_HUMAN_INHIBITOR_A/model.pdb')

# superpose all to apo
zkk = zkk.superpose(apo)
ij = ij.superpose(apo)
inhibitor = inhibitor.superpose(apo)

# cut-off appropriate parts - we're cutting off C-terminal domain PLUS 2 last residues of the SET domain
zkk_selection = zkk.top.select('not (resid 145 or resid 146 or resid 147 or resid 148 or resid 149 or resid 150 or resid 151 or resid 152 or resid 153 or resid 154 or resid 155 or resid 156 or resid 157 or resid 158 or resid 159 or resid 160 or resid 161)')
zkk_cut = zkk.atom_slice(zkk_selection)

ij_selection = ij.top.select('not (resid 145 or resid 146 or resid 147 or resid 148 or resid 149 or resid 150 or resid 151 or resid 152 or resid 153 or resid 154 or resid 155 or resid 156 or resid 157 or resid 158 or resid 159 or resid 160 or resid 161)')
ij_cut = ij.atom_slice(ij_selection)

inhibitor_selection = inhibitor.top.select('not (resid 145 or resid 146 or resid 147 or resid 148 or resid 149 or resid 150 or resid 151 or resid 152 or resid 153 or resid 154 or resid 155 or resid 156 or resid 157 or resid 158 or resid 159 or resid 160 or resid 161)')
inhibitor_cut = inhibitor.atom_slice(inhibitor_selection)

# from apo take the C-terminal domain PLUS 2 last res of SET domain
apo_selection = apo.top.select('resid 145 or resid 146 or resid 147 or resid 148 or resid 149 or resid 150 or resid 151 or resid 152 or resid 153 or resid 154 or resid 155 or resid 156 or resid 157 or resid 158 or resid 159 or resid 160 or resid 161')
apo_cut = apo.atom_slice(apo_selection)

# add final model - SET/I-SET from Inhibitor, C-terminal from 4IJ8 (i.e. State3-State3 - like SAH/H4, but 'synthetic')
inhibitor2 = inhibitor2.superpose(ij)
inhibitor2_selection = inhibitor2.top.select('not (resid 145 or resid 146 or resid 147 or resid 148 or resid 149 or resid 150 or resid 151 or resid 152 or resid 153 or resid 154 or resid 155 or resid 156 or resid 157 or resid 158 or resid 159 or resid 160 or resid 161)')
inhibitor2_cut = inhibitor2.atom_slice(inhibitor2_selection)
ij2_selection = ij.top.select('resid 145 or resid 146 or resid 147 or resid 148 or resid 149 or resid 150 or resid 151 or resid 152 or resid 153 or resid 154 or resid 155 or resid 156 or resid 157 or resid 158 or resid 159 or resid 160 or resid 161')
ij2_cut = ij.atom_slice(ij2_selection)

# stack fragments
zkk_apo = zkk_cut.stack(apo_cut)
ij_apo = ij_cut.stack(apo_cut)
inhibitor_apo = inhibitor_cut.stack(apo_cut)
inhibitor_ij = inhibitor2_cut.stack(ij2_cut)

zkk_apo.save('pdbs/1ZKK_apo.pdb')
ij_apo.save('pdbs/4IJ8_apo.pdb')
inhibitor_apo.save('pdbs/Inhibitor_apo.pdb')
inhibitor_ij.save('pdbs/Inhibitor_4IJ8.pdb')
