1) run master.sh (bash master.sh in interactive job) --> this will generate solvated structures
2) run equilibration scripts in equil_scripts individually per structure (qsub X.sh to run GPU batch jobs) - this will generate packaged (in OpenMM 7.0) structures
3) uninstall OpenMM7 (conda remove openmm), conda clean -iptls, conda install OpenMM6.3 (conda install openmm=6.3)
3) run repackage.py (python repackage.py in interactive job) --> this will repackage all structures in OpenMM 6.3
4) Ready for FAH!

COMMENTS

1) Why choose 2BQZ over 1ZKK for the peptide-SAM structure? Because of the peptide - it is longer in 2BQZ than in 1ZKK (both in terms of actually available residues, and those in SEQRES - 1ZKK had some missing, 2BQZ had all it was supposed to have) + it only had the methylated lysine (MLZ) in 2BQZ - we did not want to have to model that methyl in. 

DISADVANTAGE: the last residue of the peptide in 2BQZ has an Y, instead of I it should have (position 27 in Uniprot: P62805)

2) Sequence: we're going with the same as in Ensembler-setups 11707 and 11709 - residues 232-393, which is the most common SEQRES-claimed span in the PDB. Sequence below:

SRKSKAELQSEERKRIDELIESGKEEGMKIDLIDGKGRGVIATKQFSRGDFVVEYHGDLI
EITDAKKREALYAQDPSTGCYMYYFQYLSKTYCVDATRETNRLGRLINHSKCGNCQTKLH
DIDGVPHLILIASRDIAAGEELLYDYGDRSKASIEAHPWLKH

i.e. 2BQZ and 4IJ8 need the N-terminal S added to the AAs actually present

3) We tried hard to have an equilibration scheme that would never blow up, but it still sometimes does - depending on PDBFixer performance / hydrogen placement by openmm.app.modeller / water placement by openmm.app.modeller. If particular with_solvent structure blows up - re-run the script, try a new one. 

That's what we had to do in the final run here that was used for FAH runs. When some systems (the top 5) were already in INTERNAL, I wanted to equilibrate the rest of them and discovered that SET8_P and SET8_PH_SAH were blowing up. Replaced no_solvent and with_solvent of those from the 'repair run' - left temp_structures there the same as SET8_P shares the temp_structure with SET8_PH, and SET8_PH_SAH with SET8_P_SAH. Those ones worked ok, so I'm thinking the temp_structures were not the reason for failures - though they could have been because of the PDBFixer modelling of the SER0. 
