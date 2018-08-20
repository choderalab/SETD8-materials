# generate non-solvated structures
python generate_structures.py

# generate charges and parameters for SAM and SAH
if [ ! -d "parameters" ]; then
  mkdir parameters
fi

# cd into gaff_parameters because antechamber spits out a few files to the current dir
cd parameters

antechamber -i ../no_solvent/SAM.pdb -fi pdb -o SAM.mol2 -fo mol2 -c bcc -nc 1
parmchk2 -i SAM.mol2 -f mol2 -o frcmod.SAM

antechamber -i ../no_solvent/SAH.pdb -fi pdb -o SAH.mol2 -fo mol2 -c bcc
parmchk2 -i SAH.mol2 -f mol2 -o frcmod.SAH

# cd out into main dir
cd ..

python after_antechamber.py

# convert the ffptm.IN to ffptm.xml
python files/processAmberForceField.py files/ffptm.in files/parm99.dat files/frcmod.ff99SBildn > parameters/ffptm.xml 

python solvate.py
