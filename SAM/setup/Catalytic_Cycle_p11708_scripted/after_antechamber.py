import parmed
from io import StringIO

params_sam = parmed.amber.AmberParameterSet.from_leaprc(StringIO(u'SAM = loadMol2 parameters/SAM.mol2\nloadAmberParams parameters/frcmod.SAM'))
params_sam = parmed.openmm.OpenMMParameterSet.from_parameterset(params_sam)
params_sam.write('parameters/SAM.xml')

params_sah = parmed.amber.AmberParameterSet.from_leaprc(StringIO(u'SAH = loadMol2 parameters/SAH.mol2\nloadAmberParams parameters/frcmod.SAH'))
params_sah = parmed.openmm.OpenMMParameterSet.from_parameterset(params_sah)
params_sah.write('parameters/SAH.xml')

# convert GAFF here for reproducibility etc. rather than taking it converted from my OpenMM conversion which has not been merged yet - gaff.dat was copied into files/ from AmberTools16
# gotta set write_unused to True
params_gaff = parmed.amber.AmberParameterSet.from_leaprc(StringIO(u'loadAmberParams gaff.dat'))
params_gaff = parmed.openmm.OpenMMParameterSet.from_parameterset(params_gaff)
params_gaff.write('parameters/gaff.xml', write_unused=True)
