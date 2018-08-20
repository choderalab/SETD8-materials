# run this script to run everything!

python generate_structures.py

ce -BUILD structures/SET8_noH.pdb > structures/SET8.pdb
reduce -BUILD structures/SET8_P_noH.pdb > structures/SET8_P.pdb
reduce -BUILD structures/SET8_MeP_noH.pdb > structures/SET8_MeP.pdb
reduce -BUILD structures/SAM_noH.pdb > structures/SAM.pdb
reduce -BUILD structures/SAH_noH.pdb > structures/SAH.pdb

python after_reduce.py
