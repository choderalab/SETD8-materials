import mdtraj as md

# we are supplying the models produced in project 11707 with Ensembler from respective crystal structures 
# all are of sequence:  SRKSKAELQSEERKRIDELIESGKEEGMKIDLIDGKGRGVIATKQFSRGDFVVEYHGDLIEITDAKKREALYAQDPSTGCYMYYFQYLSKTYCVDATRETNRLGRLINHSKCGNCQTKLHDIDGVPHLILIASRDIAAGEELLYDYGDRSKASIEAHPWLKH
ij = md.load_pdb('files/SETD8_HUMAN_4IJ8_A/model.pdb')
apo = md.load_pdb('files/SETD8_HUMAN_APO_A/model.pdb')
inhibitor2 = md.load_pdb('files/KMT5A_HUMAN_TDIZ_A/model.pdb')

# superpose all to apo
ij = ij.superpose(apo)
inhibitor2 = inhibitor2.superpose(apo)

# cut-off appropriate parts - we're cutting off C-terminal domain PLUS 2 last residues of the SET domain
ij_selection = ij.top.select('not (resid 145 or resid 146 or resid 147 or resid 148 or resid 149 or resid 150 or resid 151 or resid 152 or resid 153 or resid 154 or resid 155 or resid 156 or resid 157 or resid 158 or resid 159 or resid 160 or resid 161)')
ij_cut = ij.atom_slice(ij_selection)

apo_selection = apo.top.select('not (resid 145 or resid 146 or resid 147 or resid 148 or resid 149 or resid 150 or resid 151 or resid 152 or resid 153 or resid 154 or resid 155 or resid 156 or resid 157 or resid 158 or resid 159 or resid 160 or resid 161)')
apo_cut = apo.atom_slice(apo_selection)

# from inhibitor2 take the C-terminal domain PLUS 2 last res of SET domain
inhibitor2_selection = inhibitor2.top.select('resid 145 or resid 146 or resid 147 or resid 148 or resid 149 or resid 150 or resid 151 or resid 152 or resid 153 or resid 154 or resid 155 or resid 156 or resid 157 or resid 158 or resid 159 or resid 160 or resid 161')
inhibitor2_cut = inhibitor2.atom_slice(inhibitor2_selection)

# stack fragments
ij_inhibitor2 = ij_cut.stack(inhibitor2_cut)
apo_inhibitor2 = apo_cut.stack(inhibitor2_cut)

ij_inhibitor2.save('pdbs/4IJ8_inhibitor2.pdb')
apo_inhibitor2.save('pdbs/Apo_inhibitor2.pdb')
